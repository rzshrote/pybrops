import numpy
import math
import copy

from . import CrossoverOperator

class UniformSetCrossoverOperator(CrossoverOperator):
    """docstring for UniformSetCrossoverOperator."""

    def __init__(self, M, rng, **kwargs):
        """
        M : float
        rng :
        """
        super(UniformSetCrossoverOperator, self).__init__(**kwargs)
        self.M = M
        self.rng = rng

    def nparent():
        doc = "The nparent property."
        def fget(self):
            return 2
        def fset(self, value):
            raise RuntimeError("variable 'nparent' is read-only")
        def fdel(self):
            raise RuntimeError("variable 'nparent' is read-only")
        return locals()
    nparent = property(**nparent())

    def crossover(self, parents, **kwargs):
        # copy parent data, so it is not overwritten
        ind = [copy.copy(parents[0]), copy.copy(parents[1])]

        clen = len(ind[0])                      # get chromosome length
        cdtype = ind[0].dtype                   # get chromosome dtype
        out = numpy.empty(clen, cdtype)         # create new offspring chromosome
        xoprob = numpy.empty(clen, 'float64')   # generate crossover probability array
        xoprob[0] = 0.5
        xoprob[0:] = 0.5 * (1.0 - math.exp(-2.0 * self.M / clen))

        # generate random number array
        rnd = self.rng.uniform(0, 1, clen)

        # crossover algorithm
        p = 0                                                       # get starting individual phase index
        for i in range(clen):                                       # for each point in the chromosome
            if rnd[i] < xoprob[i]:                                  # if a crossover has occured
                p = 1 - p                                           # switch parent
            while ind[p][i] in out[:i]:                             # while the candidate parent element is in offspring
                for j in range(i):                                  # for each element in offspring
                    if out[j] == ind[p][i]:                         # search if offspring element equals parent element
                        ind[p][i], ind[p][j] = ind[p][j], ind[p][i] # change element positions on parent
            out[i] = ind[p][i]                                      # save element into offspring

        return out
