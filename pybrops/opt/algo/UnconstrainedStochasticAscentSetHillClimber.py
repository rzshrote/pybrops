"""
Module implementing a stochastic ascent hill climber algorithm adapted for
subset selection optimization.
"""

import numpy

from pybrops.opt.algo.UnconstrainedOptimizationAlgorithm import UnconstrainedOptimizationAlgorithm

class UnconstrainedStochasticAscentSetHillClimber(UnconstrainedOptimizationAlgorithm):
    """
    Class implementing a stochastic ascent hill climber algorithm adapted for
    subset selection optimization. The search space is discrete and nominal in
    nature.
    """

    ########################## Special Object Methods ##########################
    def __init__(self, k, setspace, rng, objwt = 1.0, **kwargs: dict):
        """
        Constructor for a stochastic ascent set hill-climber.

        Parameters
        ----------
        k : int
            Number of states to select in set space. Select n choose k.
        setspace : numpy.ndarray
            Array of shape (n,) where 'n' is the size of the set space.
        rng : numpy.random.Generator, numpy.random.RandomState
            Random number generator.
        objwt : float, numpy.ndarray
            Weight to apply to objective function.
            If the function is a maximizing function, provide a positive weight.
            If the function is a minimizing function, provide a negative weight.
        kwargs : dict
            Additional keyword arguments.
        """
        super(UnconstrainedStochasticAscentSetHillClimber, self).__init__(**kwargs)
        self.k = k
        self.setspace = setspace
        self.rng = rng
        self.objwt = objwt

    ############################## Object Methods ##############################
    def optimize(self, objfn, k = None, setspace = None, objwt = None, **kwargs: dict):
        """
        Optimize an objective function.

        Parameters
        ----------
        objfn : callable
            Objective function to optimize.

        Returns
        -------
        out : dict
            A dictionary containing the best solution and its objective
            function evaluation.
        """
        # get parameters for optimization
        if k is None:
            k = self.k
        if setspace is None:
            setspace = self.setspace
        if objwt is None:
            objwt = self.objwt

        # initialize
        wrkss = setspace.copy()    # copy set space
        self.rng.shuffle(wrkss)    # shuffle search space

        # get starting solution and score
        gbest_soln = wrkss[:k].copy()                   # copy solution
        gbest_score = objfn(gbest_soln, **kwargs)       # score solution
        gbest_wscore = numpy.dot(gbest_score, objwt)    # weight scores

        # create exchange indices
        exchix = numpy.array(
            [[i,j] for i in range(k) for j in range(k, len(wrkss))]
        )

        # establish stopping criterion
        iterate = True

        # main loop
        while iterate:
            self.rng.shuffle(exchix)                    # shuffle exchange indices
            local_optima = True                         # whether we went through all exchange indices
            for i,j in exchix:                          # for each exchange pair
                wrkss[i], wrkss[j] = wrkss[j], wrkss[i] # exchange values
                score = objfn(wrkss[:k], **kwargs)      # score configuration
                wscore = numpy.dot(score, objwt)        # weigh score
                if wscore > gbest_wscore:               # if weighted score is better
                    gbest_soln = wrkss[:k].copy()       # copy new best solution
                    gbest_score = score                 # store new best score
                    gbest_wscore = wscore               # store new best weighted score
                    local_optima = False                # indicate we are not in a local optima
                    break                               # break from for loop
                wrkss[i], wrkss[j] = wrkss[j], wrkss[i] # exchange values back to original
            iterate = not local_optima                  # update whether to continue climbing

        # return results
        out = {
            "soln" : gbest_soln,
            "objfn_eval" : gbest_score,
            "objfn_weval" : gbest_wscore
        }

        return out
