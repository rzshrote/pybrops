import numpy

from . import OptimizationAlgorithm

class SteepestAscentSetHillClimber(OptimizationAlgorithm):
    """docstring for SteepestAscentSetHillClimber."""

    def __init__(self, k, setspace, rng, objwt = 1.0, **kwargs):
        """
        objwt : float, numpy.ndarray
            Weight to apply to objective function.
            If the function is a maximizing function, provide a positive weight.
            If the function is a minimizing function, provide a negative weight.
        """
        super(SteepestAscentSetHillClimber, self).__init__(**kwargs)
        self.k = k
        self.setspace = setspace
        self.rng = rng
        self.objwt = objwt

    def optimize(self, objfn, k = None, setspace = None, objwt = None, **kwargs):
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

        # establish stopping criterion
        iterate = True

        # main loop
        while iterate:
            exchlen = len(wrkss) - k                        # exchange length
            soln = numpy.tile(gbest_soln, (k*exchlen,1))    # tile new solutions
            for i in range(k):                              # for each dimension
                soln[i*exchlen:(i+1)*exchlen,i] = wrkss[k:] # exchange values

            # score new configurations using objfn
            score = numpy.apply_along_axis(objfn, 1, soln, **kwargs)

            # weight scores
            wscore = numpy.dot(score, objwt)

            # if any weighted scores are better than best score
            if numpy.any(wscore > gbest_score):
                gbest_ix = wscore.argmax()              # get best index
                gbest_soln = soln[gbest_ix,:].copy()    # get best solution
                gbest_score = score[gbest_ix]           # get best score
                gbest_wscore = wscore[gbest_ix]         # get best weighted score
                # reorder and update the working set space
                wrkss[:k] = gbest_soln
                mask = ~numpy.in1d(setspace, gbest_soln)    # invert
                wrkss[k:] = setspace[mask]
            else:
                iterate = False     # no better results; found a local optima

        # return results
        out = {
            "soln" : gbest_soln,
            "objfn_eval" : gbest_score,
            "objfn_weval" : gbest_wscore
        }

        return out

    def optimize_vec(fn):
        raise NotImplementedError("method is abstract")