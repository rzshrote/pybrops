"""
Module implementing an NSGA-II genetic algorithm adapted for subset selection
optimization.
"""

from typing import Union
import numpy
import math
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
from deap import benchmarks

from pybrops.core.random.prng import global_prng
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.core.util.pareto import is_pareto_efficient
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.core.error.error_type_python import check_is_int
from pybrops.core.error.error_type_python import check_is_float
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState

class NSGA2SetGeneticAlgorithm(OptimizationAlgorithm):
    """
    Class implementing an NSGA-II genetic algorithm adapted for subset selection
    optimization. The search space is discrete and nominal in nature.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, ngen = 250, mu = 100, lamb = 100, M = 1.5, rng = global_prng, **kwargs: dict):
        """
        Constructor for NSGA-II set optimization algorithm.

        Parameters
        ----------
        ngen : int
            Number of generations to evolve population.
        mu : int
            Number of parental candidates to keep in population.
        lamb : int
            Number of progeny to generate.
        M : float
            Length of the chromosome genetic map, in Morgans.
        rng : numpy.random.Generator, numpy.random.RandomState, None
            Random number generator source.
        kwargs : dict
            Additional keyword arguments.
        """
        super(NSGA2SetGeneticAlgorithm, self).__init__(**kwargs)
        self.ngen = ngen
        self.mu = mu
        self.lamb = lamb
        self.M = M
        self.rng = rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def ngen(self) -> int:
        """Number of generations."""
        return self._ngen
    @ngen.setter
    def ngen(self, value: int) -> None:
        """Set number of generations."""
        check_is_int(value, "ngen")     # must be int
        check_is_gt(value, "ngen", 0)   # int must be >0
        self._ngen = value
    @ngen.deleter
    def ngen(self) -> None:
        """Delete number of generations."""
        del self._ngen

    @property
    def mu(self) -> int:
        """Number of individuals in the main chromosome population."""
        return self._mu
    @mu.setter
    def mu(self, value: int) -> None:
        """Set number of individuals in the main chromosome population."""
        check_is_int(value, "mu")   # must be int
        check_is_gt(value, "mu", 0) # int must be >0
        self._mu = value
    @mu.deleter
    def mu(self) -> None:
        """Delete number of individuals in the main chromosome population."""
        del self._mu

    @property
    def lamb(self) -> int:
        """Number of progenies to generate from the main chromosome population."""
        return self._lamb
    @lamb.setter
    def lamb(self, value: int) -> None:
        """Set number of progenies to generate from the main chromosome population."""
        check_is_int(value, "lamb")     # must be int
        check_is_gt(value, "lamb", 0)   # int must be >0
        self._lamb = value
    @lamb.deleter
    def lamb(self) -> None:
        """Delete number of progenies to generate from the main chromosome population."""
        del self._lamb

    @property
    def M(self) -> float:
        """Length of the genetic map in Morgans."""
        return self._M
    @M.setter
    def M(self, value: float) -> None:
        """Set length of the genetic map in Morgans."""
        check_is_float(value, "M")  # must be int
        check_is_gt(value, "M", 0)  # int must be >0
        self._M = value
    @M.deleter
    def M(self) -> None:
        """Delete length of the genetic map in Morgans."""
        del self._M

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

    # define set crossover operator
    def cxSet(self, ind1, ind2, indpb):
        """
        Set crossover operator.

        Parameters
        ----------
        ind1 : numpy.ndarray
            Chromosome of the first parent (modified in place).
        ind2 : numpy.ndarray
            Chromosome of the second parent (modified in place).
        indpb : float
            Probability of initiating a crossover at a specific chromosome
            locus.

        Returns
        -------
        out : tuple
            A tuple of length 2 containing progeny resulting from the crossover.
        """
        mab = ~numpy.isin(ind1,ind2)        # get mask for ind1 not in ind2
        mba = ~numpy.isin(ind2,ind1)        # get mask for ind2 not in ind1
        ap = ind1[mab]                      # get reduced ind1 chromosome
        bp = ind2[mba]                      # get reduced ind2 chromosome
        clen = min(len(ap), len(bp))        # get minimum chromosome length
        # crossover algorithm
        p = 0                               # get starting individual phase index
        for i in range(clen):               # for each point in the chromosome
            if self.rng.random() < indpb:   # if a crossover has occured
                p = 1 - p                   # switch parent
            if p == 1:                      # if using second parent
                ap[i], bp[i] = bp[i], ap[i] # exchange alleles
        ind1[mab] = ap                      # copy over exchanges
        ind2[mba] = bp                      # copy over exchanges
        return ind1, ind2

    # define set mutation operator
    def mutSet(self, ind, sspace, indpb):
        """
        Set mutation operator.

        Parameters
        ----------
        ind : numpy.ndarray
            Individual chromosome to mutate (modified in place).
        sspace : numpy.ndarray
            Array representing the set search space.
        indpb : float
            Probability of mutation at a single locus.

        Returns
        -------
        ind : array_like
            A mutated chromosome.
        """
        for i in range(len(ind)):
            if self.rng.random() < indpb:
                mab = ~numpy.isin(sspace,ind)
                ind[i] = self.rng.choice(sspace[mab], 1, False)[0]
        return ind

    def optimize(self, objfn, k, sspace, objfn_wt, **kwargs: dict):
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

        Returns
        -------
        out : tuple
            A tuple of length 3 containing ``(frontier, sel_config, misc)``.

            Where:

            - ``frontier`` is a matrix of the Pareto frontier points.
            - ``sel_config`` is an array of corresponding decision variables
              for the corresponding Pareto frontier points.
            - ``misc`` is a dictionary of miscellaneous output.
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
        creator.create("Individual", numpy.ndarray, fitness=creator.FitnessMax)

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
        pareto_mask = is_pareto_efficient(
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
