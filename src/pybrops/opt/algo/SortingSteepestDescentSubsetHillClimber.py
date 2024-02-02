"""
Module implementing a steepest ascent hill climber capable of handling constraints.
"""

__all__ = [
    "SteepestDescentSubsetHillClimber",
]

from typing import Optional
import numpy
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.opt.algo.SubsetOptimizationAlgorithm import SubsetOptimizationAlgorithm
from pybrops.opt.prob.SubsetProblem import SubsetProblem, check_SubsetProblem_is_single_objective, check_is_SubsetProblem
from pybrops.opt.soln.SubsetSolution import SubsetSolution

class SortingSteepestDescentSubsetHillClimber(SubsetOptimizationAlgorithm):
    """
    A variant on the steepest descent hill climber for subset search spaces.
    The hillclimber first tests each element in the decision space, identifies 
    the best decision space variables, then iteratively improves the solution 
    using a steepest descent algorithm to improve the score if any contraint 
    violations exist.

    The initial sorting component of the algorithm assumes a convex search space
    where decision variables are additive.
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

    ############################ Object Properties #############################

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

        ### Phase 1: Solution by sorting
        ### ============================

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

        # get the top (lowest) ``prob.ndecn`` indices and their decision space values
        gbest_ix = ix[0:prob.ndecn,0]

        # construct an initial best solution from top elements
        gbest_soln = prob.decn_space[gbest_ix]

        ### Phase 2: Iteratively improve to reduce constraint violations
        ### ============================================================

        # get search space elements not in current global best solution
        wrkss = prob.decn_space[numpy.logical_not(numpy.in1d(prob.decn_space, gbest_soln))]

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
                # for each element not in the solution vector, but in the decision space
                for j in range(len(wrkss)):
                    # exchange values to propose new solution
                    gbest_soln[i], wrkss[j] = wrkss[j], gbest_soln[i]
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
                    gbest_soln[i], wrkss[j] = wrkss[j], gbest_soln[i]
            if (best_i is None) or (best_j is None):
                break
            gbest_soln[best_i], wrkss[best_j] = wrkss[best_j], gbest_soln[best_i] # exchange values
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