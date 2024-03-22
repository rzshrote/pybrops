"""
Module implementing an optimization algorithm that scores individual values in 
the search space, sorts the scores, and selects a solution using the sorted list.
"""

__all__ = [
    "SortingSubsetOptimizationAlgorithm",
]

from typing import Optional
import numpy
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.opt.algo.SubsetOptimizationAlgorithm import SubsetOptimizationAlgorithm
from pybrops.opt.prob.SubsetProblem import SubsetProblem
from pybrops.opt.prob.SubsetProblem import check_SubsetProblem_is_single_objective
from pybrops.opt.prob.SubsetProblem import check_is_SubsetProblem
from pybrops.opt.soln.SubsetSolution import SubsetSolution

class SortingSubsetOptimizationAlgorithm(SubsetOptimizationAlgorithm):
    """
    Optimization algorithm class that scores individuals separate from each other,
    sorts to scores, and selects a solution from the sorted list.

    Assumes a convex search space where decision variables are additive.
    Ignores any and all constraint violations.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
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
        pass

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
        prob : SubsetProblem
            A subset problem definition object on which to optimize.
        miscout : dict
            Miscellaneous output from the constrained optimizaiont algorithm.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : SubsetSolution
            An object containing the solution to the provided subset problem.
        """
        # check inputs
        check_is_SubsetProblem(prob, "prob")
        check_SubsetProblem_is_single_objective(prob, "prob")
        if miscout is not None:
            check_is_dict(miscout, "miscout")

        # evaluate each decision space element individually
        evals = [prob.evalfn(numpy.array([e])) for e in prob.decn_space]

        # unpack list of 3-tuples to a 3-tuple of Tuple[ndarray]
        obj, ineqcv, eqcv = zip(*evals)

        # convert Tuple[ndarray] to ndarray
        obj = numpy.stack(obj)          # (n,1)
        ineqcv = numpy.stack(ineqcv)    # (n,nineqcv)
        eqcv = numpy.stack(eqcv)        # (n,neqcv)
        
        # we only care about the obj; disregard the constraint violations
        # sort the objectives and get their indices
        # (n,1)
        ix = obj.argsort(0)

        # get the top (lowest) ``ndecn`` indices and their decision space values
        ndecn = prob.ndecn
        gbest_ix = ix[0:ndecn,0]
        gbest_soln = prob.decn_space[gbest_ix]

        # score the solution collectively
        gbest_obj, gbest_ineqcv, gbest_eqcv = prob.evalfn(gbest_soln)

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

        return out
