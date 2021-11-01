import numpy
import math
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
from deap import benchmarks

import pybropt.core.random

from . import OptimizationAlgorithm

class NSGA2SetGeneticAlgorithm(OptimizationAlgorithm):
    """docstring for NSGA2SetGeneticAlgorithm."""

    def __init__(self, ngen = 250, mu = 100, lamb = 100, M = 1.5, rng = None, **kwargs):
        super(NSGA2SetGeneticAlgorithm, self).__init__(**kwargs)
        self.ngen = ngen
        self.mu = mu
        self.lamb = lamb
        self.M = M
        self.rng = pybropt.core.random if rng is None else rng

    # define set crossover operator
    def cxSet(self, ind1, ind2, indpb):
        a = numpy.array(ind1)           # convert ind1 to numpy.ndarray
        b = numpy.array(ind2)           # convert ind2 to numpy.ndarray
        mab = ~numpy.isin(a,b)          # get mask for ind1 not in ind2
        mba = ~numpy.isin(b,a)          # get mask for ind2 not in ind1
        ap = a[mab]                     # get reduced ind1 chromosome
        bp = b[mba]                     # get reduced ind2 chromosome
        clen = min(len(ap), len(bp))    # get minimum chromosome length
        # crossover algorithm
        p = 0                               # get starting individual phase index
        for i in range(clen):               # for each point in the chromosome
            if self.rng.random() < indpb:   # if a crossover has occured
                p = 1 - p                   # switch parent
            if p == 1:                      # if using second parent
                ap[i], bp[i] = bp[i], ap[i] # exchange alleles
        a[mab] = ap                         # copy over exchanges
        b[mba] = bp                         # copy over exchanges
        # copy to original arrays
        for i in range(len(ind1)):
            ind1[i] = a[i]
        for i in range(len(ind2)):
            ind2[i] = b[i]
        return ind1, ind2

    # define set mutation operator
    def mutSet(self, ind, sspace, indpb):
        a = numpy.array(sspace)
        for i in range(len(ind)):
            if self.rng.random() < indpb:
                mab = ~numpy.isin(a,ind)
                ind[i] = self.rng.choice(a[mab], 1, False)[0]
        return ind

    # based on:
    # https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
    @staticmethod
    def is_pareto_efficient(fmat, wt, return_mask = True):
        """
        Find the pareto-efficient points (maximizing function)

        Parameters
        ----------
        fmat : numpy.ndarray
            A matrix of shape (npt, nobj) containing fitness values. Where
            'npt' is the number of points and 'nobj' is the number of
            objectives.
        return_mask : bool
            If True, return a mask.

        Returns
        -------
        out : numpy.ndarray
            An array of indices of pareto-efficient points.
            If return_mask is True, this will be an (npt, ) boolean array
            Otherwise it will be a (n_efficient_points, ) integer array of indices.
        """
        fmat = fmat * (wt.flatten()[None,:])    # apply weights
        npt = fmat.shape[0]                     # get number of points
        is_efficient = numpy.arange(npt)        # starting list of efficient points (holds indices)
        pt_ix = 0  # Next index in the is_efficient array to search for
        while pt_ix < len(fmat):
            ndpt_mask = numpy.any(fmat > fmat[pt_ix], axis=1)
            ndpt_mask[pt_ix] = True
            is_efficient = is_efficient[ndpt_mask]  # Remove dominated points
            fmat = fmat[ndpt_mask]
            pt_ix = numpy.sum(ndpt_mask[:pt_ix])+1
        if return_mask:
            is_efficient_mask = numpy.zeros(npt, dtype = bool)
            is_efficient_mask[is_efficient] = True
            return is_efficient_mask
        else:
            return is_efficient

    def optimize(self, objfn, k, sspace, objfn_wt, **kwargs):
        """
        Optimize an objective function.

        Parameters
        ----------
        objfn : callable
            Objective function which to optimize.
        k : int
            Number of decision variables in the search space.
            A vector is formed as sspace^k
        sspace : numpy.ndarray
            Search space that the OptimizationAlgorithm searches in.
        objfn_wt : numpy.ndarray
            Weight(s) applied to output(s) from the objfn.
        """
        # convert objective function weights to a numpy array
        if not hasattr(objfn_wt, "__iter__"):
            objfn_wt = [objfn_wt]
        objfn_wt = numpy.array(objfn_wt)

        # MOGS objectives are minimizing (DEAP uses larger fitness as better)
        creator.create(
            "FitnessMax",
            base.Fitness,
            weights = tuple(objfn_wt)
        )

        # create an individual, which is a list representation
        creator.create("Individual", list, fitness=creator.FitnessMax)

        # create a toolbox
        toolbox = base.Toolbox()

        # since this is a subset problem, represent individual as a permutation
        toolbox.register(
            "permutation",          # create a permutation protocol
            self.rng.choice,        # randomly sample
            sspace,                 # from sspace
            size = k,               # select k from sspace
            replace = False         # no replacement
        )

        # register individual creation protocol
        toolbox.register(
            "individual",           # name of function in toolbox
            tools.initIterate,
            creator.Individual,
            toolbox.permutation
        )

        # register population creation protocol
        toolbox.register(
            "population",           # name of function in toolbox
            tools.initRepeat,       # list of toolbox.individual objects
            list,                   # containter type
            toolbox.individual      # function to execute
        )

        # register the objective function
        toolbox.register(
            "evaluate",             # name of function in toolbox
            objfn                   # function to execute
        )

        ### register the crossover operator
        indpb = 0.5 * (1.0 - math.exp(-2.0 * self.M / k))
        toolbox.register(
            "mate",                 # name of function
            self.cxSet,             # function to execute
            indpb = indpb           # average of M crossovers
        )

        # register the mutation operator
        toolbox.register(
            "mutate",               # name of function
            self.mutSet,            # custom mutation operator
            sspace = sspace,        # set space
            indpb = 2.0 / k         # probability of mutation
        )

        # register the selection operator
        toolbox.register(
            "select",
            tools.selNSGA2
        )

        # register logbook statistics to take
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("min", numpy.min, axis=0)
        stats.register("max", numpy.max, axis=0)
        stats.register("avg", numpy.mean, axis=0)
        stats.register("std", numpy.std, axis=0)

        # create logbook
        logbook = tools.Logbook()
        logbook.header = ("gen", "evals", "min", "max", "avg", "std")

        # create population
        pop = toolbox.population(n = self.mu)

        # evaluate individuals with an invalid fitness
        invalid_ind = [ind for ind in pop if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # assign crowding distances (no actual selection is done)
        pop = toolbox.select(pop, len(pop))

        # compile population statistics
        record = stats.compile(pop)

        # record statistics in logbook
        logbook.record(gen=0, evals=len(invalid_ind), **record)

        # main genetic algorithm loop
        for gen in range(1, self.ngen):
            # create lambda progeny using distance crowding tournament selection
            offspring = tools.selTournamentDCD(pop, self.lamb)

            # clone individuals for modification
            offspring = [toolbox.clone(ind) for ind in offspring]

            # for each offspring pair
            for ind1, ind2 in zip(offspring[::2], offspring[1::2]):
                toolbox.mate(ind1, ind2)    # mating is guaranteed
                toolbox.mutate(ind1)        # mutate ind1
                toolbox.mutate(ind2)        # mutate ind2
                del ind1.fitness.values     # delete fitness value since ind1 is modified
                del ind2.fitness.values     # delete fitness value since ind1 is modified

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            # Select the next generation population
            pop = toolbox.select(pop + offspring, self.mu)

            # save logs
            record = stats.compile(pop)
            logbook.record(gen = gen, evals = len(invalid_ind), **record)

        # extract population fitnesses (solutions)
        pop_soln = numpy.array([ind.fitness.values for ind in pop])

        # extract population decision configurations
        pop_decn = numpy.array(pop)

        # get pareto frontier mask
        pareto_mask = self.is_pareto_efficient(
            pop_soln,
            wt = objfn_wt,
            return_mask = True
        )

        # get pareto frontier
        frontier = pop_soln[pareto_mask]

        # get selection configurations
        sel_config = pop_decn[pareto_mask]

        # stuff everything else into a dictionary
        misc = {
            "pop_soln" : pop_soln,
            "pop_decn" : pop_decn,
            "pop" : pop,
            "logbook" : logbook
        }

        return frontier, sel_config, misc
