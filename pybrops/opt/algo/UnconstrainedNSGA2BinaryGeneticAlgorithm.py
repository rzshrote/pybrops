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
from pybrops.opt.algo.UnconstrainedOptimizationAlgorithm import UnconstrainedOptimizationAlgorithm
from pybrops.core.util.pareto import is_pareto_efficient
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.core.error.error_type_python import check_is_int
from pybrops.core.error.error_type_python import check_is_float
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState

class UnconstrainedNSGA2BinaryGeneticAlgorithm(UnconstrainedOptimizationAlgorithm):
    """
    Class implementing an NSGA-II genetic algorithm adapted for subset selection
    optimization. The search space is binary in nature.
    """

    ########################## Special Object Methods ##########################
    def __init__(self, ngen = 250, mu = 100, lamb = 100, rng = global_prng, **kwargs: dict):
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
        super(UnconstrainedNSGA2BinaryGeneticAlgorithm, self).__init__(**kwargs)
        self.ngen = ngen
        self.mu = mu
        self.lamb = lamb
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

    ############################## Object Methods ##############################

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

        # since this is a subset problem, represent individual as a binary decision
        toolbox.register(
            "binary",               # create a binary protocol
            self.rng.choice,        # randomly sample
            [0,1],                  # from sspace
            size = k,               # select k from sspace
            replace = True          # with replacement
        )

        # register individual creation protocol
        toolbox.register(
            "individual",           # name of function in toolbox
            tools.initIterate,
            creator.Individual,
            toolbox.binary
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
        toolbox.register(
            "mate",                 # name of function
            tools.cxTwoPoint        # function to execute
        )

        # register the mutation operator
        toolbox.register(
            "mutate",               # name of function
            tools.mutFlipBit,       # bit flip mutation
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
