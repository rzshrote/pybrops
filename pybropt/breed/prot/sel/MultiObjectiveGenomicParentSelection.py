import numpy
import math
import types

from deap import algorithms
from deap import base
from deap import creator
from deap import tools
from deap import benchmarks

import pybropt.core.random

from . import SelectionProtocol

from pybropt.core.error import check_is_int
from pybropt.core.error import check_is_str

class MultiObjectiveGenomicParentSelection(SelectionProtocol):
    """docstring for MultiObjectiveGenomicParentSelection."""

    def __init__(self, nparent, ncross, nprogeny, algorithm, method,
    objfn_trans, objfn_trans_kwargs = None, objfn_wt = 1.0,
    ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = 1.0,
    target = "positive", weight = "magnitude",
    ga_ngen = 250, ga_mu = 100, ga_lamb = 100, ga_M = 1.5,
    rng = None, **kwargs):
        """
        nparent : int
            Number of parents to select.
        ncross : int
            Number of crosses per configuration.
        nprogeny : int
            Number of progeny to derive from each cross.
        algorithm : Algorithm
            Optimization algorithm to optimize the objective function.
        method : str
            Method of selecting parents.
            Method   | Description
            ---------+----------------------------------------------------------
            "single" | MOGS is transformed to a single objective and
                     | optimization is done on the transformed function.
                     | This is done using the 'trans' function provided:
                     |    optimize : objfn_trans(MOGS)
            -------+------------------------------------------------------------
            "pareto" | MOGS is transformed by a transformation function, but NOT
                     | reduced to a single objective. The Pareto frontier for
                     | this transformed function is mapped using a
                     | multi-objective GA.
                     | Objectives are scaled to [0,1] and a vector orthogonal to
                     | the hyperplane defined by the extremes of the front is
                     | drawn starting at the point defined by 'ndset_trans'. The
                     | closest point on the Pareto frontier to the orthogonal
                     | vector is selected.
            ---------+----------------------------------------------------------
        objfn_trans : function or callable
            Function to transform the MOGS function. If method = "single", this
            function must return a scalar. If method = "pareto", this function
            must return a numpy.ndarray.

            Function definition:
            --------------------
            objfn_trans(obj, **kwargs):
                Parameters
                    obj : scalar, numpy.ndarray
                        Objective scalar or vector to be transformed
                    **kwargs : dict
                        Additional keyword arguments
                Returns
                    out : scalar, numpy.ndarray
                        Transformed objective scalar or vector.
        objfn_trans_kwargs : dict
            Dictionary of keyword arguments to be passed to 'objfn_trans'.
        objfn_wt : float, numpy.ndarray
            Weight applied to transformed objective function. Indicates whether
            a function is maximizing or minimizing.
                1.0 for maximizing function.
                -1.0 for minimizing function.
        ndset_trans : numpy.ndarray
            Function to transform nondominated points along the Pareto frontier
            into a single score for each point.

            Function definition:
            --------------------
            ndset_trans(ndset, **kwargs):
                Parameters
                    ndset : numpy.ndarray
                        Array of shape (j, o) containing nondominated points.
                        Where 'j' is the number of nondominated points and 'o'
                        is the number of objectives.
                    **kwargs : dict
                        Additional keyword arguments.
                Returns
                    out : numpy.ndarray
                        Array of shape (j,) containing transformed Pareto
                        frontier points.
        ndset_trans_kwargs : dict
            Dictionary of keyword arguments to be passed to 'ndset_trans'.
        ndset_wt : float
            Weight applied to transformed nondominated points along Pareto
            frontier. Indicates whether a function is maximizing or minimizing.
                1.0 for maximizing function.
                -1.0 for minimizing function.
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
        ga_ngen : int
            Number of generations to evolve multi-objective genetic algorithm.
        ga_mu : int
            Number of individuals in the main population for the multi-objective
            genetic algorithm.
        ga_lamb : int
            Number of progenies to generate per generation for the multi-
            objective genetic algorithm.
        ga_M : float
            Length of genetic map to simulate crossover. Units are in Morgans.
            Number of recombinations follows a Poisson distribution meaning that
            the mean number of recombinations per chromosome in the genetic
            algorithm is 'ga_M'.
            ONLY RELEVANT TO THE GENETIC ALGORITHM. DOES NOT AFFECT SIMULATION
            OF GENETIC RECOMBINATIONS IN BREEDING SIMULATIONS.
        rng : numpy.random.Generator or None
            A random number generator source. Used for optimization algorithms.
            If 'rng' is None, use pybropt.core.random module (NOT THREAD SAFE!).
        """
        super(MultiObjectiveGenomicParentSelection, self).__init__(**kwargs)

        # error checks
        check_is_int(nparent, "nparent")
        check_is_int(ncross, "ncross")
        check_is_int(nprogeny, "nprogeny")
        # check_is_Algorithm(algorithm, "algorithm") # TODO: implement me
        check_is_str(method, "method")
        check_is_int(ga_ngen, "ga_ngen")
        check_is_int(ga_mu, "ga_mu")
        check_is_int(ga_lamb, "ga_lamb")

        # variable assignment
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.algorithm = algorithm
        self.method = method
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = {} if objfn_trans_kwargs is None else objfn_trans_kwargs
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = {} if ndset_trans_kwargs is None else ndset_trans_kwargs
        self.ndset_wt = ndset_wt
        self.target = target
        self.weight = weight
        self.ga_ngen = ga_ngen
        self.ga_mu = ga_mu
        self.ga_lamb = ga_lamb
        self.ga_M = ga_M
        self.rng = pybropt.core.random if rng is None else rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def pselect(self, t_cur, t_max, geno, bval, gmod,
    nparent = None, algorithm = None, method = None,
    objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None,
    ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = None,
    **kwargs):
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
        nparent : int
            Number of parents to select
        **kwargs
            Additional keyword arguments to be passed to either
            algorithm.optimize (method = "single") or self.ppareto (method =
            "pareto"), depending on method used.

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
        # if any optional parameters are None, set to defaults.
        if nparent is None:
            nparent = self.nparent
        if algorithm is None:
            algorithm = self.algorithm
        if method is None:
            method = self.method
        if objfn_trans is None:
            objfn_trans = self.objfn_trans
        if objfn_trans_kwargs is None:
            objfn_trans_kwargs = self.objfn_trans_kwargs
        if objfn_wt is None:
            objfn_wt = self.objfn_wt
        if ndset_trans is None:
            ndset_trans = self.ndset_trans
        if ndset_trans_kwargs is None:
            ndset_trans_kwargs = self.ndset_trans_kwargs
        if ndset_wt is None:
            ndset_wt = self.ndset_wt

        # convert method string to lower
        method = method.lower()

        # protocols for which method to use
        if method == "single":
            # get objective function
            objfn = self.pobjfn(
                t_cur = t_cur,
                t_max = t_max,
                geno = geno,
                bval = bval,
                gmod = gmod,
                trans = objfn_trans,
                trans_kwargs = objfn_trans_kwargs
            )

            # optimize the objective function
            soln_dict = algorithm.optimize(
                objfn,
                k = nparent,
                setspace = numpy.arange(geno["cand"].ntaxa),
                objwt = objfn_wt
            )

            # extract solution
            sel = soln_dict["soln"]

            # return optimized selection configuration and other items
            return geno["cand"], sel, self.ncross, self.nprogeny, soln_dict

        elif method == "pareto":
            # get the pareto frontier
            frontier, sel_config, misc = self.ppareto(
                t_cur = t_cur,
                t_max = t_max,
                geno = geno,
                bval = bval,
                gmod = gmod,
                nparent = nparent,
                objfn_trans = objfn_trans,
                objfn_trans_kwargs = objfn_trans_kwargs,
                objfn_wt = objfn_wt
            )

            # get scores for each of the points along the pareto frontier
            score = ndset_wt * ndset_trans(frontier, **ndset_trans_kwargs)

            # get index of maximum score
            ix = score.argmax()

            # add fields to misc
            misc["frontier"] = frontier
            misc["sel_config"] = sel_config

            return geno["cand"], sel_config[ix], self.ncross, self.nprogeny, misc
        else:
            raise ValueError("argument 'method' must be either 'single' or 'pareto'")

    def pobjfn(self, t_cur, t_max, geno, bval, gmod, trans, trans_kwargs, **kwargs):
        """
        Return a parent selection objective function.

        Parameters
        ----------
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:
                Must accept a single numpy.ndarray argument.
                Must return a single object, whether scalar or numpy.ndarray.
        **kwargs : dict
            Dictionary of keyword arguments to pass to 'trans' function.
        """
        mat = geno["cand"].mat                      # (m,n,p) get genotype matrix
        beta = gmod["cand"].beta                    # (p,t) get regression coefficients
        mkrwt = self.calc_mkrwt(self.weight, beta)  # (p,t) get marker weights
        tfreq = self.calc_tfreq(self.target, beta)  # (p,t) get target allele frequencies

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn.__code__,                        # byte code pointer
            self.objfn.__globals__,                     # global variables
            None,                                       # new name for the function
            (mat, tfreq, mkrwt, trans, trans_kwargs),   # default values for arguments
            self.objfn.__closure__                      # closure byte code pointer
        )

        return outfn

    # TODO: implementation of this function
    def pobjfn_vec(self, t_cur, t_max, geno, bval, gmod, **kwargs):
        """
        Return a vectorized objective function.
        """
        raise NotImplementedError("method is abstract")

    def ppareto(self, t_cur, t_max, geno, bval, gmod, nparent = None,
    objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None,
    ngen = None, mu = None, lamb = None, M = None, **kwargs):
        # get parameters
        if nparent is None:
            nparent = self.nparent
        if objfn_trans is None:
            objfn_trans = self.objfn_trans
        if objfn_trans_kwargs is None:
            objfn_trans_kwargs = self.objfn_trans_kwargs
        if objfn_wt is None:
            objfn_wt = self.objfn_wt
        if ngen is None:
            ngen = self.ga_ngen
        if mu is None:
            mu = self.ga_mu
        if lamb is None:
            lamb = self.ga_lamb
        if M is None:
            M = self.ga_M

        # get number of taxa
        ntaxa = geno["cand"].ntaxa

        # get objective function
        objfn = self.pobjfn(
            t_cur = t_cur,
            t_max = t_max,
            geno = geno,
            bval = bval,
            gmod = gmod,
            trans = objfn_trans,
            trans_kwargs = objfn_trans_kwargs,
            **kwargs
        )

        # MOGS objectives are minimizing (DEAP uses larger fitness as better)
        creator.create("FitnessMax", base.Fitness, weights = tuple(objfn_wt))

        # create an individual, which is a list representation
        creator.create("Individual", list, fitness=creator.FitnessMax)

        # create a toolbox
        toolbox = base.Toolbox()

        # since this is a subset problem, represent individual as a permutation
        toolbox.register(
            "permutation",          # create a permutation protocol
            self.rng.choice,        # randomly sample
            ntaxa,                  # from 0:ntaxa
            size = nparent,         # select nparent from 0:ntaxa
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
        indpb = 0.5 * (1.0 - math.exp(-2.0 * M / nparent))
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
            indpb = 2.0 / nparent                   # probability of mutation
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
        pop = toolbox.population(n = mu)

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

        # extract population fitnesses (solutions)
        pop_solution = numpy.array([ind.fitness.values for ind in pop])

        # extract population decision configurations
        pop_decision = numpy.array(pop)

        # get pareto frontier mask
        pareto_mask = self.is_pareto_efficient(
            pop_solution,
            wt = objfn_wt,
            return_mask = True
        )

        # get pareto frontier
        frontier = pop_solution[pareto_mask]

        # get selection configurations
        sel_config = pop_decision[pareto_mask]

        # stuff everything else into a dictionary
        misc = {
            "pop_solution" : pop_solution,
            "pop_decision" : pop_decision,
            "pop" : pop,
            "logbook" : logbook
        }

        return frontier, sel_config, misc

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn(sel, mat, tfreq, mkrwt, trans, kwargs):
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
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:
                Must accept a single numpy.ndarray argument.
                Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to 'trans' function.

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

        # transform and return
        return trans(mogs, **kwargs)

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

    @staticmethod
    def calc_mkrwt(weight, beta):
        if isinstance(weight, str):
            weight = weight.lower()             # convert to lowercase
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

    @staticmethod
    def calc_tfreq(target, beta):
        if isinstance(target, str):
            target = target.lower()                 # convert to lowercase
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

    @staticmethod
    def traitsum_trans(mat, wt = None):
        if wt is None:
            return mat.sum(-1)
        else:
            return mat @ wt

    @staticmethod
    def weightsum_trans(mat, wt = None):
        if wt is None:
            return mat.sum()
        else:
            return mat.dot(wt)

    @staticmethod
    def vecptdist_trans(mat, objfn_wt, wt = None):
        """
        mat : numpy.ndarray
            (npt, nobj)
        objfn_wt : numpy.ndarray
            (nobj,)
        wt : numpy.ndarray
            (nobj,)
        """
        # create a default wt if wt is None
        if wt is None:
            wt = numpy.ones(mat.shape[1], dtype = 'float64')

        # transform mat to all maximizing functions
        # (nobj,) -> (1,nobj)
        # (npt,nobj) * (1,nobj) -> (npt,nobj)
        mat = mat * wt[None,:]

        # subtract column minimums
        mat = mat - mat.min(0)

        # divide by column maximums; mat is in range [0,1]
        mat = mat / mat.max(0)

        # calculate distance between point and line
        # calculate ((v dot p) / (v dot v)) * v
        # where v is the line vector originating from 0
        #       p is the point vector

        # get inverse of (v dot v)
        # (nobj,) dot (nobj,) -> scalar
        vdvinv = 1.0 / objfn_wt.dot(objfn_wt)

        # get scaling factor (v dot p) / (v dot v)
        # (npt,nobj) dot (nobj,) -> (npt,)
        # (npt,) * scalar -> (npt,)
        scale = mat.dot(objfn_wt) * vdvinv

        # use outer product to get points on plane that intersect line
        # (npt,) outer (nobj,) -> (npt,nobj)
        P = numpy.outer(scale, objfn_wt)

        # calculate difference between each vector
        # (npt,nobj) - (npt,nobj) -> (npt,nobj)
        diff = mat - P

        # take vector norms of difference to get distances
        # (npt,nobj) -> (npt,)
        d = numpy.linalg.norm(diff, axis = 1)

        return d
