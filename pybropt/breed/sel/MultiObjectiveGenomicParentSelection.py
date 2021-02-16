import numpy
import math

from deap import algorithms
from deap import base
from deap import creator
from deap import tools
from deap import benchmarks


from . import ParentSelectionOperator

from pybropt.core.error import check_is_int

class MultiObjectiveGenomicParentSelection(ParentSelectionOperator):
    """docstring for MultiObjectiveGenomicParentSelection."""

    def __init__(self, k_p, traitobjwt_p, traitsum_p, objsum_p, algorithm_p, ncross, nprogeny, rng, target = "positive", weight = "magnitude", **kwargs):
        """
        k_p : int
            Number of parents to select.
        target : str or numpy.ndarray
            If target is a string, check value and follow these rules:
                Value         | Description
                --------------+-------------------------------------------------
                "positive"    | Select alleles with the most positive effect.
                "negative"    | Select alleles with the most negate effect.
                "stabilizing" | Set target allele frequency to 0.5.
            If target is a numpy.ndarray, use values as is.
        weight : str or numpy.ndarray
            If weight is a string, check value and follow these rules:
                Value       | Description
                ------------+---------------------------------------------------
                "magnitude" | Assign weights using the magnitudes of regression coefficients.
                "equal"     | Assign weights equally.
        """
        super(MultiObjectiveGenomicParentSelection, self).__init__(**kwargs)

        # error checks
        check_is_int(k_p, "k_p")

        # variable assignment
        self.k_p = k_p
        self.traitobjwt_p = traitobjwt_p
        self.traitsum_p = traitsum_p
        self.objsum_p = objsum_p
        self.algorithm_p = algorithm_p
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.rng = rng
        self.target = target
        self.weight = weight

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def calc_mkrwt(self, weight, beta):
        if isinstance(weight, str):
            if weight == "magnitude":           # return abs(beta)
                return numpy.absolute(beta)
            elif weight == "equal":             # return 1s matrix
                return numpy.full(beta.shape, 1.0, dtype='float64')
            else:
                raise ValueError("string value for 'weight' not recognized")
        elif isinstance(weight, numpy.ndarray):
            return weight
        else:
            raise TypeError("variable 'weight' must be a string or numpy.ndarray")

    def calc_tfreq(self, target, beta):
        if isinstance(target, str):
            if target == "positive":
                return numpy.float64(beta >= 0.0)   # positive alleles are desired
            elif target == "negative":
                return numpy.float64(beta <= 0.0)   # negative alleles are desired
            elif target == "stabilizing":
                return 0.5                          # both alleles desired
                # return numpy.full(coeff.shape, 0.5, dtype = 'float64')
            else:
                raise ValueError("string value for 'target' not recognized")
        elif isinstance(target, numpy.ndarray):
            return target
        else:
            raise TypeError("variable 'target' must be a string or numpy.ndarray")

    def pselect(self, t_cur, t_max, geno, bval, gmod, k = None, traitobjwt = None, traitsum = None, objsum = None, algorithm = None, **kwargs):
        """
        Select parents individuals for breeding.

        Parameters
        ----------
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        geno : dict
            A dict containing genotypic data for all breeding populations.
            Must have the following fields:
                Field | Type                         | Description
                ------+------------------------------+--------------------------
                cand  | PhasedGenotypeMatrix         | Parental candidate breeding population
                main  | PhasedGenotypeMatrix         | Main breeding population
                queue | List of PhasedGenotypeMatrix | Breeding populations on queue
                ""
        bval : dict
            A dict containing breeding value data.
            Must have the following fields:
                Field      | Type                        | Description
                -----------+-----------------------------+----------------------
                cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                main       | BreedingValueMatrix         | Main breeding population breeding values
                main_true  | BreedingValueMatrix         | Main breeding population true breeding values
        gmod : dict
            A dict containing genomic models.
            Must have the following fields:
                Field | Type                 | Description
                ------+----------------------+----------------------------------
                cand  | GenomicModel         | Parental candidate breeding population genomic model
                main  | GenomicModel         | Main breeding population genomic model
                true  | GenomicModel         | True genomic model for trait(s)
        k : int
        traitobjwt : None, numpy.ndarray
            Combined trait weights and objective weights matrix.
            This matrix must be compatible with shape (2,t).
            If None, default to settings established at object construction.
            Specifying:
                Only trait weights: (1,t)
                Only objective weights: (2,1)
                Trait and objective weights: (2,t)
        traitsum : None, bool
            Sum across traits.
            If None, default to settings established at object construction.
            If True:
                (2,t) -> (2,)
            Else:
                (2,t)
        objsum : None, bool
            Sum across objectives.
            If None, default to settings established at object construction.
            If True:
                (2,t) -> (t,)
            Else:
                (2,t)
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing five objects: (pgvmat, sel, ncross, nprogeny, misc)
            pgvmat : PhasedGenotypeVariantMatrix
                A PhasedGenotypeVariantMatrix of parental candidates.
            sel : numpy.ndarray
                Array of indices specifying a cross pattern. Each index
                corresponds to an individual in 'pgvmat'.
            ncross : numpy.ndarray
                Number of crosses to perform per cross pattern.
            nprogeny : numpy.ndarray
                Number of progeny to generate per cross.
            misc : dict
                Miscellaneous output (user defined).
        """
        # get parameters
        if k is None:
            k = self.k_p
        if traitobjwt is None:
            traitobjwt = self.traitobjwt_p
        if traitsum is None:
            traitsum = self.traitsum_p
        if objsum is None:
            objsum = self.objsum_p
        if algorithm is None:
            algorithm = self.algorithm_p

        # get objective function
        objfn = self.pobjfn(
            t_cur = t_cur,
            t_max = t_max,
            geno = geno,
            bval = bval,
            gmod = gmod,
            traitobjwt = traitobjwt,
            traitsum = traitsum,
            objsum = objsum
        )

        soln_dict = algorithm.optimize(
            objfn,
            k = k,
            setspace = numpy.arange(geno["cand"].ntaxa),
            objwt = 1.0
        )

        # extract solution
        sel = soln_dict["soln"]

        misc = soln_dict

        return geno["cand"], sel, self.ncross, self.nprogeny, misc

    def pobjfn(self, t_cur, t_max, geno, bval, gmod, traitobjwt, traitsum, objsum, **kwargs):
        """
        Return a parent selection objective function.

        Parameters
        ----------
        traitobjwt : numpy.ndarray, None
            Combined trait weights and objective weights matrix.
            This matrix must be compatible with shape (2,t).
            Specifying:
                Only trait weights: (1,t)
                Only objective weights: (2,1)
                Trait and objective weights: (2,t)
        traitsum : bool
            Sum across traits.
            If True:
                (2,t) -> (2,)
            Else:
                (2,t)
        objsum : bool
            Sum across objectives.
            If True:
                (2,t) -> (t,)
            Else:
                (2,t)
        """
        mat = geno["cand"].mat                      # (m,n,p) get genotype matrix
        beta = gmod["cand"].beta                    # (p,t) get regression coefficients
        mkrwt = self.calc_mkrwt(self.weight, beta)  # (p,t) get marker weights
        tfreq = self.calc_tfreq(self.target, beta)  # (p,t) get target allele frequencies

        def objfn(sel, mat = mat, tfreq = tfreq, mkrwt = mkrwt, traitobjwt = traitobjwt, traitsum = traitsum, objsum = objsum):
            """
            Multi-objective genomic selection objective function.
                The goal is to minimize this function. Lower is better.
                This is a bare bones function. Minimal error checking is done.

            Given a 2D weight vector 'dcoeff', calculate the Euclidian distance from the
            origin according to:
                dist = dot( dcoeff, F(x) )
                Where:
                    F(x) is a vector of objective functions:
                        F(x) = < f_PAU(x), f_PAFD(x) >

            f_PAU(x):

            Given the provided genotype matrix 'geno' and row selections from it 'sel',
            calculate the selection allele freq. From the selection allele frequencies
            and the target allele frequencies, determine if the target frequencies
            cannot be attained after unlimited generations and selection rounds.
            Multiply this vector by a weight coefficients vector 'wcoeff'.

            f_PAFD(x):

            Given a genotype matrix, a target allele frequency vector, and a vector of
            weights, calculate the distance between the selection frequency and the
            target frequency.

            Parameters
            ==========
            sel : numpy.ndarray, None
                A selection indices matrix of shape (k,)
                Where:
                    'k' is the number of individuals to select.
                Each index indicates which individuals to select.
                Each index in 'sel' represents a single individual's row.
                If 'sel' is None, use all individuals.
            mat : numpy.ndarray, None
                A int8 binary genotype matrix of shape (m, n, p).
                Where:
                    'm' is the number of chromosome phases (2 for diploid, etc.).
                    'n' is the number of individuals.
                    'p' is the number of markers.
                Remarks:
                    Shape of the matrix is most critical. Underlying matrix
                    operations will support other numeric data types.
            tfreq : floating, numpy.ndarray
                A target allele frequency matrix of shape (p, t).
                Where:
                    'p' is the number of markers.
                    't' is the number of traits.
                Example:
                    tfreq = numpy.array([0.2, 0.6, 0.7])
            mkrwt : numpy.ndarray
                A marker weight coefficients matrix of shape (p, t).
                Where:
                    'p' is the number of markers.
                    't' is the number of traits.
                Remarks: Values in 'mkrwt' have an assumption:
                    All values must be non-negative.
            traitobjwt : numpy.ndarray, None
                Combined trait weights and objective weights matrix.
                This matrix must be compatible with shape (2,t).
                Specifying:
                    Only trait weights: (1,t)
                    Only objective weights: (2,1)
                    Trait and objective weights: (2,t)
            traitsum : bool
                Sum across traits.
                If True:
                    (2,t) -> (2,)
                Else:
                    (2,t)
            objsum : bool
                Sum across objectives.
                If True:
                    (2,t) -> (t,)
                Else:
                    (2,t)

            Returns
            =======
            mogs : numpy.ndarray
                A MOGS score matrix of shape (2,t) or other.
            """
            # if no selection, select all
            if sel is None:
                sel = slice(None)

            # generate a view of the genotype matrix that only contains 'sel' rows.
            # (m,(k,),p) -> (m,k,p)
            sgeno = mat[:,sel,:]

            # calculate reciprocal number of phases
            rphase = 1.0 / (sgeno.shape[0] * sgeno.shape[1])

            # calculate population frequencies; add axis for correct broadcast
            # (m,k,p).sum[0,1] -> (p,)
            # (p,) * scalar -> (p,)
            # (p,None) -> (p,1) # need (p,1) for broadcasting with (p,t) arrays
            pfreq = (sgeno.sum((0,1)) * rphase)[:,None]

            # calculate some inequalities for use multiple times
            pfreq_lteq_0 = (pfreq <= 0.0)   # is population freq <= 0.0
            pfreq_gteq_1 = (pfreq >= 1.0)   # is population freq >= 1.0

            # calculate allele unavailability
            allele_unavail = numpy.where(
                tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
                pfreq_lteq_0,           # then set True if sel has allele freq == 0
                numpy.where(            # else
                    tfreq > 0.0,        # if 0.0 < target freq < 1.0
                    numpy.logical_or(   # then set True if pop freq is outside (0.0,1.0)
                        pfreq_lteq_0,
                        pfreq_gteq_1
                    ),
                    pfreq_gteq_1        # else set True if pop freq is >= 1.0
                )
            )

            # calculate distance between target and population
            # (p,t)-(p,1) -> (p,t)
            dist = numpy.absolute(tfreq - pfreq)

            # compute f_PAU(x)
            # (p,t) * (p,t) -> (p,t)
            # (p,t).sum[0] -> (t,)
            pau = (mkrwt * allele_unavail).sum(0)

            # compute f_PAFD(x)
            # (p,t) * (p,t) -> (p,t)
            # (p,t).sum[0] -> (t,)
            pafd = (mkrwt * dist).sum(0)

            # stack to make MOGS matrix
            # (2,t)
            mogs = numpy.stack([pau, pafd])

            # weight traits
            if traitobjwt is not None:
                mogs *= traitobjwt

            # sum across any axes as needed.
            if objsum or traitsum:
                axistuple = tuple()
                if objsum:
                    axistuple += (0,)
                if traitsum:
                    axistuple += (1,)
                mogs = mogs.sum(axistuple)

            return mogs

        return objfn

    # TODO: implementation of this function
    # def pobjfn_vec(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
    #     """
    #     Return a vectorized objective function.
    #     """
    #     mat = geno["cand"].mat      # genotype matrix
    #     mu = gmod["cand"].mu        # trait means
    #     beta = gmod["cand"].beta    # regression coefficients
    #
    #     def objfn_vec(sel, mat = mat, mu = mu, beta = beta, traitwt = traitwt):
    #         """
    #         Score a population of individuals based on Conventional Genomic Selection
    #         (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
    #         Genomic Estimated Breeding Values (GEBV) for a population.
    #
    #         Parameters
    #         ----------
    #         sel : numpy.ndarray
    #             A selection indices matrix of shape (j,k)
    #             Where:
    #                 'j' is the number of selection configurations.
    #                 'k' is the number of individuals to select.
    #             Each index indicates which individuals to select.
    #             Each index in 'sel' represents a single individual's row.
    #             If 'sel' is None, use all individuals.
    #         mat : numpy.ndarray
    #             A int8 binary genotype matrix of shape (m, n, p).
    #             Where:
    #                 'm' is the number of chromosome phases (2 for diploid, etc.).
    #                 'n' is the number of individuals.
    #                 'p' is the number of markers.
    #         mu : numpy.ndarray
    #             A trait mean matrix of shape (t, 1)
    #             Where:
    #                 't' is the number of traits.
    #         beta : numpy.ndarray
    #             A trait prediction coefficients matrix of shape (p, t).
    #             Where:
    #                 'p' is the number of markers.
    #                 't' is the number of traits.
    #         traitwt : numpy.ndarray, None
    #             A trait objective coefficients matrix of shape (t,).
    #             Where:
    #                 't' is the number of objectives.
    #             These are used to weigh objectives in the weight sum method.
    #             If None, do not multiply GEBVs by a weight sum vector.
    #
    #         Returns
    #         -------
    #         cgs : numpy.ndarray
    #             A trait GEBV matrix of shape (j,k,t) if objwt is None.
    #             A trait GEBV matrix of shape (j,k) if objwt shape is (t,)
    #             OR
    #             A weighted GEBV matrix of shape (t,).
    #             Where:
    #                 'k' is the number of individuals selected.
    #                 't' is the number of traits.
    #         """
    #         # (m,n,p)[:,(j,k),:] -> (m,j,k,p)
    #         # (m,j,k,p) -> (j,k,p)
    #         # (j,k,p) . (p,t) -> (j,k,t)
    #         # (j,k,t) + (1,t) -> (j,k,t)
    #         cgs = mat[:,sel,:].sum(0).dot(beta) + mu.T
    #
    #         # (j,k,t) . (t,) -> (j,k)
    #         if traitwt is not None:
    #             cgs = cgs.dot(traitwt)
    #
    #         return cgs
    #
    #     return objfn_vec

    # https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
    def is_pareto_efficient(self, fitnesses, return_mask = True):
        """
        Find the pareto-efficient points (maximizing function)
        :param fitnesses: An (n_points, n_fitnesses) array
        :param return_mask: True to return a mask
        :return: An array of indices of pareto-efficient points.
            If return_mask is True, this will be an (n_points, ) boolean array
            Otherwise it will be a (n_efficient_points, ) integer array of indices.
        """
        is_efficient = numpy.arange(fitnesses.shape[0])
        n_points = fitnesses.shape[0]
        next_point_index = 0  # Next index in the is_efficient array to search for
        while next_point_index<len(fitnesses):
            nondominated_point_mask = numpy.any(fitnesses>fitnesses[next_point_index], axis=1)
            nondominated_point_mask[next_point_index] = True
            is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
            fitnesses = fitnesses[nondominated_point_mask]
            next_point_index = numpy.sum(nondominated_point_mask[:next_point_index])+1
        if return_mask:
            is_efficient_mask = numpy.zeros(n_points, dtype = bool)
            is_efficient_mask[is_efficient] = True
            return is_efficient_mask
        else:
            return is_efficient

    def ppareto(self, t_cur, t_max, geno, bval, gmod, k = None, traitobjwt = None, traitsum = None, ngen = 250, mu = 100, lamb = 100, M = 1.5, **kwargs):
        # get parameters
        if k is None:
            k = self.k_p
        if traitobjwt is None:
            traitobjwt = self.traitobjwt_p
        if traitsum is None:
            traitsum = self.traitsum_p

        # get number of taxa
        ntaxa = geno["cand"].ntaxa

        # get objective function
        objfn = self.pobjfn(
            t_cur = t_cur,
            t_max = t_max,
            geno = geno,
            bval = bval,
            gmod = gmod,
            traitobjwt = traitobjwt,
            traitsum = traitsum,
            objsum = False
        )

        # MOGS objectives are minimizing (DEAP uses larger fitness as better)
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,-1.0))

        # create an individual, which is a list representation
        creator.create("Individual", list, fitness=creator.FitnessMin)

        # create a toolbox
        toolbox = base.Toolbox()

        # since this is a subset problem, represent individual as a permutation
        toolbox.register(
            "permutation",          # create a permutation protocol
            self.rng.choice,        # randomly sample
            ntaxa,                  # from 0:ntaxa
            size = k,               # select k from 0:ntaxa
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

        # define set crossover operator
        def cxSet(ind1, ind2, indpb):
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

        ### register the crossover operator
        indpb = 0.5 * (1.0 - math.exp(-2.0 * M / k))
        toolbox.register(
            "mate",                 # name of function
            cxSet,                  # function to execute
            indpb = indpb           # average of M crossovers
        )

        # define set mutation operator
        def mutSet(ind, setspace, indpb):
            a = numpy.array(setspace)
            for i in range(len(ind)):
                if self.rng.random() < indpb:
                    mab = ~numpy.isin(a,ind)
                    ind[i] = self.rng.choice(a[mab], 1, False)[0]
            return ind

        # register the mutation operator
        toolbox.register(
            "mutate",                               # name of function
            mutSet,                                 # custom mutation operator
            setspace = [i for i in range(ntaxa)],   # set space
            indpb = 2.0 / k                         # probability of mutation
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
        logbook.header = "gen", "evals", "min", "max", "avg", "std"

        # create population
        pop = toolbox.population(n = mu)

        # evaluate individuals with an invalid fitness
        invalid_ind = [ind for ind in pop if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # assign crowding distances (no actual selection is done)
        pop = toolbox.select(pop, len(pop))

        # for i,ind in enumerate(pop):
        #     print("ind", i, ":", ind)

        # compile population statistics
        record = stats.compile(pop)

        # record statistics in logbook
        logbook.record(gen=0, evals=len(invalid_ind), **record)

        # main genetic algorithm loop
        for gen in range(1, ngen):
            # create lambda progeny using distance crowding tournament selection
            offspring = tools.selTournamentDCD(pop, lamb)

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
            pop = toolbox.select(pop + offspring, mu)

            # save logs
            record = stats.compile(pop)
            logbook.record(gen = gen, evals = len(invalid_ind), **record)

        # extract frontier points
        frontier = numpy.array([ind.fitness.values for ind in pop])

        return frontier, pop, logbook
