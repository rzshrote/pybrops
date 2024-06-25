"""
Module implementing a steepest ascent hill climber capable of handling constraints.
"""

__all__ = [
    "NeighborSteepestDescentSubsetHillClimber",
]

from numbers import Integral
from typing import Callable, Optional
from typing import Union
import numpy
from numpy.random import Generator
from numpy.random import RandomState
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_type_python import check_is_Callable, check_is_dict, check_is_str
from pybrops.core.random.prng import global_prng
from pybrops.opt.algo.SubsetOptimizationAlgorithm import SubsetOptimizationAlgorithm
from pybrops.opt.prob.SubsetProblem import SubsetProblem
from pybrops.opt.prob.SubsetProblem import check_SubsetProblem_is_single_objective
from pybrops.opt.prob.SubsetProblem import check_is_SubsetProblem
from pybrops.opt.soln.SubsetSolution import SubsetSolution

def matrix_coordinate_neighbors(index: Integral, grid: numpy.ndarray, exclude: numpy.ndarray):
    # get the grid coordinates for the element
    # (s,d) -> (d,)
    ecoords = grid[index]
    # for each dimension, get indices of neighbors, concatenate
    nidx = numpy.concatenate([numpy.flatnonzero(grid[:,i] == ecoords[i]) for i in range(len(ecoords))])
    # remove any from exclude list
    nidx = nidx[~numpy.in1d(nidx, exclude)]
    return nidx

def matrix_coordinate_neighbors_with_random(index: Integral, grid: numpy.ndarray, exclude: numpy.ndarray, nrand: int = 100):
    # get the grid coordinates for the element
    # (s,d) -> (d,)
    ecoords = grid[index]
    # construct neighborhood indices, eliminate excluded indices
    neighbors_idx = numpy.unique(numpy.concatenate([numpy.flatnonzero(grid[:,i] == ecoords[i]) for i in range(len(ecoords))]))
    neighbors_idx = neighbors_idx[~numpy.in1d(neighbors_idx, exclude)]
    # construct indices to exclude from random selection
    exclude_idx = numpy.unique(numpy.concatenate([neighbors_idx, exclude]))
    # construct random candidates indices
    randcand_idx = numpy.delete(numpy.arange(len(grid)), exclude_idx)
    # select random indices
    rand_idx = numpy.random.choice(randcand_idx, min(nrand, len(randcand_idx)), replace = False)
    # concatenate neighbors + random indices
    out_idx = numpy.concatenate([neighbors_idx, rand_idx])
    return out_idx

class NeighborSteepestDescentSubsetHillClimber(
        SubsetOptimizationAlgorithm,
    ):
    """
    Steepest descent hill climber for subset search spaces.
    """
    ########################## Special Object Methods ##########################
    def __init__(
            self,
            neighborattr: str,
            neighborfn: Callable,
            rng: Optional[Union[Generator,RandomState]] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for a steepest ascent hillclimber capable of handling
        constraints.
        
        Parameters
        ----------
        rng : 
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        self.neighborattr = neighborattr
        self.neighborfn = neighborfn
        self.rng = rng
    ############################ Object Properties #############################
    @property
    def neighborattr(self) -> str:
        """neighborattr."""
        return self._neighborattr
    @neighborattr.setter
    def neighborattr(self, value: str) -> None:
        """Set neighborattr."""
        check_is_str(value, "neighborattr")
        self._neighborattr = value
    @property
    def neighborfn(self) -> Callable:
        """neighborfn."""
        return self._neighborfn
    @neighborfn.setter
    def neighborfn(self, value: Callable) -> None:
        """Set neighborfn."""
        check_is_Callable(value, "neighborfn")
        self._neighborfn = value
    @property
    def rng(self) -> Union[Generator,RandomState]:
        """Random number generator source."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[Generator,RandomState]) -> None:
        """Set random number generator source."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng")
        self._rng = value
    ############################## Object Methods ##############################
    def minimize(
            self, 
            prob: SubsetProblem,
            miscout: Optional[dict] = None,
            **kwargs: dict
        ) -> SubsetSolution:
        """
        Minimize an optimization problem.

        Parameters
        ----------
        prob : SetProblem
            A problem definition object on which to optimize.
        miscout : dict
            Miscellaneous output from the constrained optimizaiont algorithm.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : SubsetSolution
            An object containing the solution to the provided problem.
        """
        # check inputs
        check_is_SubsetProblem(prob, "prob")
        check_SubsetProblem_is_single_objective(prob, "prob")
        if miscout is not None:
            check_is_dict(miscout, "miscout")
        # get grid
        # (s,d)
        grid = getattr(prob, self.neighborattr)
        # randomly initialize a solution
        gbest_soln = self.rng.choice(prob.decn_space, prob.ndecn)
        # get starting solution score and constraint violations
        gbest_obj, gbest_ineqcv, gbest_eqcv = prob.evalfn(gbest_soln)
        gbest_score = gbest_obj.sum()
        gbest_cv = gbest_ineqcv.sum() + gbest_eqcv.sum()
        # hillclimber
        while True:
            best_i = None
            best_j = None
            best_obj, best_ineqcv, best_eqcv = gbest_obj, gbest_ineqcv, gbest_eqcv
            best_score = gbest_score
            best_cv = gbest_cv
            # for each element in the solution vector
            for i in range(len(gbest_soln)):
                # get neighbors
                neighbors = self.neighborfn(gbest_soln[i], grid, gbest_soln)
                # for each neighbor not in the solution vector, but in the decision space
                for j in neighbors:
                    # exchange values to propose new solution
                    old = gbest_soln[i]
                    gbest_soln[i] = j
                    # score proposed solution
                    prop_obj, prop_ineqcv, prop_eqcv = prob.evalfn(gbest_soln)
                    prop_score = prop_obj.sum()
                    prop_cv = prop_ineqcv.sum() + prop_eqcv.sum()
                    # determine if the proposed solution is better
                    # always prefer less constraint violation first
                    if prop_cv < best_cv:
                        best_i = i
                        best_j = j
                        best_obj, best_ineqcv, best_eqcv = prop_obj, prop_ineqcv, prop_eqcv
                        best_score = prop_score
                        best_cv = prop_cv
                    # if constraint violations are identical, prefer better (min) score
                    elif (prop_cv == best_cv) and (prop_score < best_score):
                        best_i = i
                        best_j = j
                        best_obj, best_ineqcv, best_eqcv = prop_obj, prop_ineqcv, prop_eqcv
                        best_score = prop_score
                        best_cv = prop_cv
                    # exchange values back to original solution
                    gbest_soln[i] = old
            if (best_i is None) or (best_j is None):
                break
            # exchange best value
            gbest_soln[best_i] = best_j
            gbest_obj, gbest_ineqcv, gbest_eqcv = best_obj, best_ineqcv, best_eqcv
            gbest_score = best_score
            gbest_cv = best_cv
        # create output
        out = SubsetSolution(
            ndecn = prob.ndecn,
            decn_space = prob.decn_space,
            decn_space_lower = prob.decn_space_lower,
            decn_space_upper = prob.decn_space_upper,
            nobj = prob.nobj,
            obj_wt = prob.obj_wt,
            nineqcv = prob.nineqcv,
            ineqcv_wt = prob.ineqcv_wt,
            neqcv = prob.neqcv,
            eqcv_wt = prob.eqcv_wt,
            nsoln = 1,
            soln_decn = numpy.stack([gbest_soln]),
            soln_obj = numpy.stack([gbest_obj]),
            soln_ineqcv = numpy.stack([gbest_ineqcv]),
            soln_eqcv = numpy.stack([gbest_eqcv]),
        )
        # add miscellaneous outputs to miscout dictionary
        if miscout is not None:
            miscout["gbest_score"] = gbest_score
            miscout["gbest_cv"] = gbest_cv
        return out

