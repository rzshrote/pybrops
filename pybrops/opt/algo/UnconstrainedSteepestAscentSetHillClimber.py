"""
Module implementing a steepest ascent hill climber algorithm adapted for subset
selection optimization.
"""

from typing import Union
import numpy

from pybrops.opt.algo.UnconstrainedOptimizationAlgorithm import UnconstrainedOptimizationAlgorithm
from pybrops.core.random.prng import global_prng
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState

class UnconstrainedSteepestAscentSetHillClimber(UnconstrainedOptimizationAlgorithm):
    """
    Class implementing a steepest ascent hill climber algorithm adapted for
    subset selection optimization. The search space is discrete and nominal in
    nature.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            rng = global_prng, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for a steepest ascent set hill-climber.

        Parameters
        ----------
        rng : numpy.random.Generator, numpy.random.RandomState
            Random number generator.
        kwargs : dict
            Additional keyword arguments.
        """
        super(UnconstrainedSteepestAscentSetHillClimber, self).__init__(**kwargs)
        self.rng = rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def rng(self) -> Union[numpy.random.Generator,numpy.random.RandomState]:
        """Random number generator source."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[numpy.random.Generator,numpy.random.RandomState]) -> None:
        """Set random number generator source."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng")
        self._rng = value
    @rng.deleter
    def rng(self) -> None:
        """Delete random number generator source."""
        del self._rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def optimize(self, objfn, k, sspace, objfn_wt, **kwargs: dict):
        # randomly initialize a solution
        gbest_soln = self.rng.choice(sspace, (k,))
        wrkss = sspace[numpy.logical_not(numpy.in1d(sspace, gbest_soln))]     # get search space elements not in chromosome

        # get starting solution score
        gbest_score = objfn(gbest_soln)
        gbest_score = gbest_score if hasattr(gbest_score, "__iter__") else (gbest_score,)
        objfn_wt = objfn_wt if hasattr(objfn_wt, "__iter__") else (objfn_wt,)
        gbest_wscore = numpy.dot(gbest_score, objfn_wt)

        # hillclimber
        while True:
            best_i = None
            best_j = None
            best_score = gbest_soln
            best_wscore = gbest_wscore
            for i in range(len(gbest_soln)):
                for j in range(len(wrkss)):
                    gbest_soln[i], wrkss[j] = wrkss[j], gbest_soln[i] # exchange values
                    score = objfn(gbest_soln)
                    score = score if hasattr(score, "__iter__") else (score,)
                    wscore = numpy.dot(score, objfn_wt)
                    if wscore > best_wscore:
                        best_i = i
                        best_j = j
                        best_score = score
                        best_wscore = wscore
                    gbest_soln[i], wrkss[j] = wrkss[j], gbest_soln[i] # exchange values back to original
            if (best_i is None) or (best_j is None):
                break
            gbest_soln[best_i], wrkss[best_j] = wrkss[best_j], gbest_soln[best_i] # exchange values
            gbest_score = best_score
            gbest_wscore = best_wscore

        # return results
        out = {
            "soln" : gbest_soln,
            "objfn_eval" : gbest_score,
            "objfn_weval" : gbest_wscore
        }

        return gbest_score, gbest_soln, out

    def optimize_vec(fn):
        raise NotImplementedError("method is abstract")
