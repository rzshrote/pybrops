import numpy

from . import MemeticOperator

class StochasticSetHillClimberMemeticOperator(MemeticOperator):
    """docstring for StochasticSetHillClimberMemeticOperator."""

    def __init__(self, setspace, rng, objwt = 1.0, **kwargs: dict):
        super(StochasticSetHillClimberMemeticOperator, self).__init__(**kwargs)
        self.setspace = setspace
        self.rng = rng
        self.objwt = objwt

    def evolve(self, objfn, pop, score, setspace = None, objwt = None, **kwargs: dict):
        """
        Parameters
        ----------
        objfn : function
        pop : list of numpy.ndarray
        score : numpy.ndarray
        kwargs : dict
            Additional arguments to pass to objfn.

        Returns
        -------
        new_pop : list of numpy.ndarray
        new_score : numpy.ndarray
        """
        if setspace is None:
            setspace = self.setspace
        if objwt is None:
            objwt = self.objwt

        new_pop = [None]*len(pop)       # list to store new individuals
        new_score = numpy.empty(        # array to store new scores
            score.shape,
            dtype = score.dtype
        )

        for c,(indiv,iscore) in enumerate(zip(pop,score)):      # index, individual, individual score
            wrkss = setspace.copy()                             # copy set space
            k = len(indiv)                                      # get length of individual
            l = len(wrkss)                                      # get length of set space
            for i in range(k):                                  # for each element in indiv
                for j in range(i,l):                            # for each element in wrkss
                    if indiv[i] == wrkss[j]:                    # if indiv element == wrkss element
                        wrkss[i], wrkss[j] = wrkss[j], wrkss[i] # exchange wrkss elements
                        break                                   # break from inner loop

            # get starting solution and score
            best_soln = indiv.copy()                            # copy starting solution array
            best_score = score                                  # get starting solution score
            best_wscore = numpy.dot(best_score, objwt)          # weight scores

            # create exchange indices
            exchix = numpy.array([[i,j] for i in range(k) for k in range(k, l)])

            # hill climber algorithm
            iterate = True                                  # stopping criterion variable
            while iterate:                                  # main hill climbing loop
                self.rng.shuffle(exchix)                    # shuffle exchange indices
                local_optima = True                         # whether we went through all exchange indices
                for i,j in exchix:                          # for each exchange pair
                    wrkss[i], wrkss[j] = wrkss[j], wrkss[i] # exchange values
                    score = objfn(wrkss[:k], **kwargs)      # score configuration
                    wscore = numpy.dot(score, objwt)        # weigh score
                    if wscore > best_wscore:                # if weighted score is better
                        best_soln = wrkss[:k].copy()        # copy new best solution
                        best_score = score                  # store new best score
                        best_wscore = wscore                # store new best weighted score
                        local_optima = False                # indicate we are not in a local optima
                        break                               # break from for loop
                    wrkss[i], wrkss[j] = wrkss[j], wrkss[i] # exchange values back to original
                iterate = not local_optima                  # update whether to continue climbing

            # store values
            new_pop[c] = best_soln
            new_score[c] = best_score

        return new_pop, new_score
