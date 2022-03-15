from math import factorial
from deap import tools

from pybrops.algo.opt.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_float
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.random import global_prng

class NSGA3UnityConstraintGeneticAlgorithm(OptimizationAlgorithm):
    """
    Class implementing an NSGA-II genetic algorithm adapted for continuous
    search space optimization with a unity constraint in the decision space.
    The search space is continuous and constrained in nature.

    A unity constraint means that genes in a chromosome representation must have
    ranges of [0.0, 1.0] and the sum of the chromosome genes must be 1.0.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, ngen = 250, mu = 100, lamb = 100, cxeta = 30.0, muteta = 20.0,
    refpnts = None, save_logbook = False, rng = None, **kwargs):
        """
        Constructor for NSGA-III optimization algorithm.

        Parameters
        ----------
        ngen : int
            Number of generations to evolve population.
        mu : int
            Number of parental candidates to keep in population.
        lamb : int
            Number of progeny to generate.
        cxeta : float
            Bounded simulated binary crossover variance parameter.
        muteta : float
            Bounded polynomial mutation variance parameter.
        refpnts : numpy.ndarray, None
            Reference points used to guide the NSGA-III algorithm through space.
        save_logbook : bool
            Whether to record a logbook during algorithm evolution and output
            the results after algorithm completion.
        rng : numpy.random.Generator, numpy.random.RandomState, None
            Random number generator source.
        kwargs : dict
            Additional keyword arguments.
        """
        super(NSGA3UnityConstraintGeneticAlgorithm, self).__init__(**kwargs)

        # perform error checks
        check_is_int(ngen, "ngen")
        check_is_int(mu, "mu")
        check_is_int(lamb, "lamb")
        check_is_float(cxeta, "cxeta")
        check_is_float(muteta, "muteta")
        if refpnts is not None:
            check_is_ndarray(refpnts, "refpnts")
        if rng is None:
            rng = global_prng
        check_is_Generator_or_RandomState(rng, "rng")

        # assign variables
        self.ngen = ngen
        self.mu = mu
        self.lamb = lamb
        self.cxeta = cxeta
        self.muteta = muteta
        self.refpnts = refpnts
        self.rng = rng

    ############################################################################
    ########################## Private Object Methods ##########################
    ############################################################################
    @staticmethod
    def _calc_H(k, x):
        """
        Calculate the number of uniform reference points the space
        dimensionality (k) and a density (x).

        Parameters
        ----------
        k : int
            Number of decision variables in the search space.
        x : int
            Density of the uniform reference points
        """
        return factorial(k+x-1) // (factorial(k-1) * factorial(x))

    def _get_refpnts(self, k):
        """
        Parameters
        ----------
        k : int
            Number of decision variables in the search space.
        """
        refpnts = self.refpnts                              # get reference points stored within self
        if refpnts is None:                                 # if no reference points provided
            H0 = self.mu - 4                                # permit 4 individuals not assigned to a reference point
            x = 1                                           # start with 1 per dimension
            while _calc_H(k, x) < H0:                       # iterate through x until we exceed H0
                x += 1                                      # increment x
            refpnts = tools.uniform_reference_points(k, x)  # create uniform reference points
        elif refpnts.shape[1] != k:                         # make sure our reference point dimensionality is compatible
            raise ValueError(
                "provided reference points are of incorrect dimensionality for problem: {0} != {1}".format(refpnts.shape[1], k)
            )
        return refpnts

    ############################################################################
    ############################## Helper Methods ##############################
    ############################################################################

    # define individual generator operator
    def unifUnitySample(self, low, up):
        """
        Uniformly sample individuals from a search space, then perform unity
        repair on the individual.

        Parameters
        ----------
        low : numpy.ndarray
            The lower bound of the search space.
        up : numpy.ndarray
            The upper bound of the search space.

        Returns
        -------
        out : list
            An individual represented as a list.
        """
        # sample from uniform search space
        out = list(self.rng.uniform(low, up))

        ### perform unity sum repair on the individual
        # calculate the inverse sum
        invsum = 1.0 / sum(out)

        # loop through and adjust chromosomes to sum to one
        for i in range(len(out)):
            out[i] *= invsum

        # return individual wrapped in tuple
        return out

    # define crossover operator
    def cxUnitySimulatedBinaryBounded(self, ind1, ind2, eta, low, up):
        """
        Perform a simulated binary crossover that modify in-place the input
        individuals, then apply a unity sum repair on the individual.

        Parameters
        ----------
        ind1 : Sequence
            The first individual participating in the crossover.
        ind2 : Sequence
            The second individual participating in the crossover.
        eta : float
            Crowding degree of the crossover. A high eta will produce children
            resembling to their parents, while a small eta will produce
            solutions much more different.
        low : float, Sequence
            The lower bound of the search space.
        up : float, Sequence
            The upper bound of the search space.

        Returns
        -------
        out : tuple
            A tuple of two individuals.
        """
        # perform bounded simulated binary crossover
        out1, out2 = tools.cxSimulatedBinaryBounded(ind1, ind2, eta, low, up)

        ### perform unity sum repair on the individual
        # calculate the inverse sum for first individual
        invsum = 1.0 / sum(out1)

        # loop through and adjust chromosomes to sum to one for first individual
        for i in range(len(out1)):
            out1[i] *= invsum

        # calculate the inverse sum for second individual
        invsum = 1.0 / sum(out2)

        # loop through and adjust chromosomes to sum to one for second individual
        for i in range(len(out2)):
            out2[i] *= invsum

        # return individuals wrapped in tuple
        return out1, out2

    # define mutation operator
    def mutUnityPolynomialBounded(self, ind, eta, low, up, indpb):
        """
        Perform a bounded polynomial mutation as implemented in the original
        NSGA-II algorithm, then apply a unity sum repair on the individual.

        Parameters
        ----------
        individual : Sequence
            Individual to be mutated.
        eta : float
            Crowding degree of the mutation. A high eta will produce a mutant
            resembling its parent, while a small eta will produce a solution
            much more different.
        low : float, Sequence
            The lower bound of the search space.
        up : float, Sequence
            The upper bound of the search space.
        indpb : float
            Probability that a gene along the chromosome will be mutated.

        Returns
        -------
        out : Sequence
            A mutated individual.
        """
        # perform bounded polynomial mutation
        out, = tools.mutPolynomialBounded(ind, eta, low, up, indpb)

        ### perform unity sum repair on the individual
        # calculate the inverse sum
        invsum = 1.0 / sum(out)

        # loop through and adjust chromosomes to sum to one
        for i in range(len(out)):
            out[i] *= invsum

        # return individual wrapped in tuple
        return out

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def optimize(self, objfn, k, sspace, objwt, **kwargs):
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
            A search space constraint array of shape ``(2,k)``. The first row
            represents the lower bound to search. The second row represents the
            upper bound to search. The first row must have elements all equal to
            ``0.0``, and the second row must have elements all equal to ``1.0``.
        objwt : numpy.ndarray
            Weight(s) applied to output(s) from the objfn.

        Returns
        -------
        out : tuple
            A tuple of length 3 (soln, decn, misc)
        """
        # convert objective function weights to a numpy array
        if not hasattr(objfn_wt, "__iter__"):
            objfn_wt = [objfn_wt]
        objfn_wt = numpy.array(objfn_wt)

        # OCS objectives are maximizing (DEAP uses larger fitness as better)
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
            "attr_float",           # create a floating point representation
            self.unifUnitySample,   # random uniform sampling
            low = sspace[0],        # search space lower bounds
            up = sspace[1]          # search space upper bounds
        )

        # register individual creation protocol
        toolbox.register(
            "individual",           # name of function in toolbox
            tools.initIterate,
            creator.Individual,
            toolbox.attr_float      # floating point chromosome generator
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
            "mate",                             # name of function
            self.cxUnitySimulatedBinaryBounded, # crossover function to execute
            eta = self.cxeta,                   # variance parameter
            low = sspace[0],                    # search space lower bounds
            up = sspace[1]                      # search space upper bounds
        )

        # register the mutation operator
        toolbox.register(
            "mutate",                       # name of function
            self.mutUnityPolynomialBounded, # custom mutation operator
            eta = self.cxeta,               # variance parameter
            low = sspace[0],                # search space lower bounds
            up = sspace[1],                 # search space upper bounds
            indpb = 1.0 / k                 # probability of mutation
        )

        # register the selection operator
        toolbox.register(
            "select",                           # name of the function
            tools.selNSGA3,                     # use NSGA-III selection operator
            ref_points = self._get_refpnts(k)   # reference point grid
        )

        # register logbook statistics to take
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("min", numpy.min, axis=0)
        stats.register("max", numpy.max, axis=0)
        stats.register("avg", numpy.mean, axis=0)
        stats.register("std", numpy.std, axis=0)

        # create logbook
        logbook = None
        if save_logbook:
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

        if save_logbook:                # optionally save logbook
            record = stats.compile(pop) # compile population statistics
            logbook.record(             # record statistics in logbook
                gen = 0,
                evals = len(invalid_ind),
                **record
            )

        # main genetic algorithm loop
        for gen in range(1, self.ngen):
            # create lambda progeny using random parent selection
            offspring = tools.selRandom(pop, self.lamb)

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
            if save_logbook:                # optionally save logbook
                record = stats.compile(pop) # compile population statistics
                logbook.record(             # record statistics in logbook
                    gen = gen,
                    evals = len(invalid_ind),
                    **record
                )

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
