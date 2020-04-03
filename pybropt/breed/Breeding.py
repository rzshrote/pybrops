# 3rd party
import numpy
import pandas

# our libraries
import pybropt.util

class Breeding:
    """
    Base class to represent a breeding strategy.
    """

    ############################################################################
    ######################### Reserved Object Methods ##########################
    ############################################################################
    def __init__(self, population, cross, method = None):
        """
        Constructor for Breeding class.

        Parameters
        ----------
        population : Population
            Breeding population object.
        cross : Cross
            Cross structure object.
        method : str, None
            Breeding method.
        """
        # initialize history lists, etc.
        self.reset()

        # type checking
        pybropt.util.check_is_Population(population, "population")
        pybropt.util.check_is_Cross(cross, "cross")
        pybropt.util.cond_check_is_string(method, "method")

        # set private variables
        self._population = population
        self._cross = cross
        self._method = method

    ############################################################################
    ################################ Properties ################################
    ############################################################################
    def population():
        doc = "The population property."
        def fget(self):
            return self._population
        def fset(self, value):
            self._population = value
        def fdel(self):
            del self._population
        return locals()
    population = property(**population())

    def cross():
        doc = "The cross property."
        def fget(self):
            return self._cross
        def fset(self, value):
            self._cross = value
        def fdel(self):
            del self._cross
        return locals()
    cross = property(**cross())

    def method():
        doc = "The method property."
        def fget(self):
            return self._method
        def fset(self, value):
            self._method = value
        def fdel(self):
            del self._method
        return locals()
    method = property(**method())

    ########################################################
    # population history properties
    def population_method():
        doc = "The population_method property."
        def fget(self):
            return self._population_method
        def fset(self, value):
            self._population_method = value
        def fdel(self):
            del self._population_method
        return locals()
    population_method = property(**population_method())

    def population_algorithm():
        doc = "The population_algorithm property."
        def fget(self):
            return self._population_algorithm
        def fset(self, value):
            self._population_algorithm = value
        def fdel(self):
            del self._population_algorithm
        return locals()
    population_algorithm = property(**population_algorithm())

    def population_seed():
        doc = "The population_seed property."
        def fget(self):
            return self._population_seed
        def fset(self, value):
            self._population_seed = value
        def fdel(self):
            del self._population_seed
        return locals()
    population_seed = property(**population_seed())

    def population_cycle():
        doc = "The population_cycle property."
        def fget(self):
            return self._population_cycle
        def fset(self, value):
            self._population_cycle = value
        def fdel(self):
            del self._population_cycle
        return locals()
    population_cycle = property(**population_cycle())

    def population_score():
        doc = "The population_score property."
        def fget(self):
            return self._population_score
        def fset(self, value):
            self._population_score = value
        def fdel(self):
            del self._population_score
        return locals()
    population_score = property(**population_score())

    def population_gebv():
        doc = "The population_gebv property."
        def fget(self):
            return self._population_gebv
        def fset(self, value):
            self._population_gebv = value
        def fdel(self):
            del self._population_gebv
        return locals()
    population_gebv = property(**population_gebv())

    def population_genotype():
        doc = "The population_genotype property."
        def fget(self):
            return self._population_genotype
        def fset(self, value):
            self._population_genotype = value
        def fdel(self):
            del self._population_genotype
        return locals()
    population_genotype = property(**population_genotype())

    ########################################################
    # selection history properties
    def selection_method():
        doc = "The selection_method property."
        def fget(self):
            return self._selection_method
        def fset(self, value):
            self._selection_method = value
        def fdel(self):
            del self._selection_method
        return locals()
    selection_method = property(**selection_method())

    def selection_algorithm():
        doc = "The selection_algorithm property."
        def fget(self):
            return self._selection_algorithm
        def fset(self, value):
            self._selection_algorithm = value
        def fdel(self):
            del self._selection_algorithm
        return locals()
    selection_algorithm = property(**selection_algorithm())

    def selection_seed():
        doc = "The selection_seed property."
        def fget(self):
            return self._selection_seed
        def fset(self, value):
            self._selection_seed = value
        def fdel(self):
            del self._selection_seed
        return locals()
    selection_seed = property(**selection_seed())

    def selection_cycle():
        doc = "The selection_cycle property."
        def fget(self):
            return self._selection_cycle
        def fset(self, value):
            self._selection_cycle = value
        def fdel(self):
            del self._selection_cycle
        return locals()
    selection_cycle = property(**selection_cycle())

    def selection_score():
        doc = "The selection_score property."
        def fget(self):
            return self._selection_score
        def fset(self, value):
            self._selection_score = value
        def fdel(self):
            del self._selection_score
        return locals()
    selection_score = property(**selection_score())

    def selection_gebv():
        doc = "The selection_gebv property."
        def fget(self):
            return self._selection_gebv
        def fset(self, value):
            self._selection_gebv = value
        def fdel(self):
            del self._selection_gebv
        return locals()
    selection_gebv = property(**selection_gebv())

    def selection_genotype():
        doc = "The selection_genotype property."
        def fget(self):
            return self._selection_genotype
        def fset(self, value):
            self._selection_genotype = value
        def fdel(self):
            del self._selection_genotype
        return locals()
    selection_genotype = property(**selection_genotype())

    ############################################################################
    ########################## Abstract Class Methods ##########################
    ############################################################################
    def objfn(self, sel, **kwargs):
        """
        Breeding method objective function. Implement this in derived classes.

        Parameters
        ----------
        sel : numpy.ndarray
            Vector of selection data.
        **kwargs : dict
            Additional arguments for this function.

        Returns
        -------
        error : NotImplementedError
            This function raises an error on execution.
        """
        raise NotImplementedError("The method 'objfn' is abstract.")

    def objfn_vec(self, sel, **kwargs):
        """
        Breeding method objective function. Implement this in derived classes.

        Parameters
        ----------
        sel : numpy.ndarray
            A 2D matrix of selection data.
        **kwargs : dict
            Additional arguments for this function.

        Returns
        -------
        error : NotImplementedError
            This function raises an error on execution.
        """
        raise NotImplementedError("The method 'objfn_vec' is abstract.")

    ############################################################################
    ########################## Private Class Methods ###########################
    ############################################################################
    def _coerce_gebv(self, gebv):
        if not isinstance(gebv, numpy.ndarray):
            # force to 2d matrix and return
            return numpy.array([[gebv]])
        elif gebv.ndim == 1:
            # add axis and return
            return numpy.array(gebv[None])
        elif gebv.ndim == 2:
            # copy and return
            return numpy.array(gebv)
        else:
            raise ValueError("'gebv' is not the correct shape")

    def _coerce_1d_string_(self, var, varname, nentry):
        if not isinstance(var, numpy.ndarray):
            # repeat and return
            return numpy.repeat(numpy.string_(str(var)), nentry)
        elif var.ndim == 1:
            if var.size == 1:
                # repeat and return
                return numpy.repeat(numpy.string_(var), nentry)
            elif var.size == nentry:
                # copy and return
                return numpy.string_([str(e) for e in var])
        raise ValueError("'%s' is not the correct shape" % varname)

    def _coerce_1d_array(self, var, varname, nentry):
        if not isinstance(var, numpy.ndarray):
            # repeat and return
            return numpy.repeat([var], nentry)
        elif var.ndim == 1:
            if var.size == 1:
                # repeat and return
                return numpy.repeat(var, nentry)
            elif var.size == nentry:
                # copy and return
                return numpy.array(var)
        raise ValueError("'%s' is not the correct shape" % varname)

    def _coerce_2d_array(self, var, varname, nentry):
        if not isinstance(var, numpy.ndarray):
            # repeat, reshape, and return
            return numpy.repeat(var, nentry).reshape(nentry, 1)
        elif var.ndim == 1:
            # tile and return
            return numpy.tile(var, (nentry, 1))
        elif var.ndim == 2 and len(var) == nentry:
            # copy and return
            return numpy.array(var)
        raise ValueError("'%s' is not the correct shape" % varname)

    def _coerce_geno(self, geno, nentry):
        if geno.ndim == 3:
            # check for correct number of genotypes
            if geno.shape[1] == nentry:
                # copy and return
                return numpy.array(geno)
        elif geno.ndim == 2:
            # reshape, copy, return
            return numpy.array(geno[:,None])
        raise ValueError("'geno' is not the correct shape")

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    ########################################################
    # population history methods
    def population_history_reset(self):
        self._population_method = []
        self._population_algorithm = []
        self._population_seed = []
        self._population_cycle = []
        self._population_score = []
        self._population_genotype = []
        self._population_gebv = []

    def add_population_history(self, method, algorithm, seed, cycle, score, gebv, geno = None):
        """
        method : str
        algorithm : str
        seed : int
        cycle : int
        score :
        """
        # must start with coercing GEBVs to 2d array
        gebv = self._coerce_gebv(gebv)

        # coerce the rest to thier correct shape
        method = self._coerce_1d_string_(method, "method", len(gebv))
        algorithm = self._coerce_1d_string_(algorithm, "algorithm", len(gebv))
        seed = self._coerce_1d_array(seed, "seed", len(gebv))
        cycle = self._coerce_1d_array(cycle, "cycle", len(gebv))
        score = self._coerce_2d_array(score, "score", len(gebv))
        if geno is not None:
            geno = self._coerce_geno(geno, len(gebv))

        self._population_method.append(method)
        self._population_algorithm.append(algorithm)
        self._population_seed.append(seed)
        self._population_cycle.append(cycle)
        self._population_score.append(score)
        self._population_gebv.append(gebv)
        if geno is not None:
            self._population_genotype.append(geno)

    def concatenate_population(self):
        if len(self._population_method) > 1:
            self._population_method = [numpy.concatenate(self._population_method, axis=0)]
        if len(self._population_algorithm) > 1:
            self._population_algorithm = [numpy.concatenate(self._population_algorithm, axis=0)]
        if len(self._population_seed) > 1:
            self._population_seed = [numpy.concatenate(self._population_seed, axis=0)]
        if len(self._population_cycle) > 1:
            self._population_cycle = [numpy.concatenate(self._population_cycle, axis=0)]
        if len(self._population_score) > 1:
            self._population_score = [numpy.concatenate(self._population_score, axis=0)]
        if len(self._population_gebv) > 1:
            self._population_gebv = [numpy.concatenate(self._population_gebv, axis=0)]
        if len(self._population_genotype) > 1:
            self._population_genotype = [numpy.concatenate(self._population_genotype, axis=1)]

    def population_history_to_dict(self, zfill = 6):
        self.concatenate_population()

        # create population stats dict
        population_dict = {
            "method" : self._population_method[0],
            "algorithm" : self._population_algorithm[0],
            "seed" : self._population_seed[0],
            "cycle" : self._population_cycle[0]
        }

        # add scores column by column
        population_dict.update(dict(zip(
            # column names
            ["score"+str(e+1) for e in range(self._population_score[0].shape[1])],
            # columns
            self._population_score[0].T
        )))

        # add GEBVs column by column
        population_dict.update(dict(zip(
            # column names
            ["gebv"+str(e+1) for e in range(self._population_gebv[0].shape[1])],
            # iterate through columns
            self._population_gebv[0].T
        )))

        # add genotypes if we have them
        if len(self._population_genotype) > 0:
            shape = self._population_genotype[0].shape

            # add first phase keys to population_dict
            population_dict.update(dict(zip(
                # column names
                ["p1_"+str(e).zfill(zfill) for e in range(shape[1])],
                # iterate through columns
                self._population_gebv[0][0].T
            )))

            # add second phase keys to population_dict
            population_dict.update(dict(zip(
                # column names
                ["p2_"+str(e).zfill(zfill) for e in range(shape[1])],
                # iterate through columns
                self._population_gebv[0][1].T
            )))

        return population_dict

    def population_history_to_df(self, zfill = 6):
        # make dictionary
        population_dict = self.population_history_to_dict(zfill)

        # make pandas.DataFrame from dict
        population_df = pandas.DataFrame(population_dict)

        return population_df

    def population_history_to_csv(self, fname, zfill = 6, **kwargs):
        population_df = self.population_history_to_df(zfill)
        population_df.to_csv(fname, **kwargs)

    ########################################################
    # selection history methods
    def selection_history_reset(self):
        self._selection_method = []
        self._selection_algorithm = []
        self._selection_seed = []
        self._selection_cycle = []
        self._selection_score = []
        self._selection_genotype = []
        self._selection_gebv = []

    def add_selection_history(self, method, algorithm, seed, cycle, score, gebv, geno = None):
        # must start with coercing GEBVs to 2d array
        gebv = self._coerce_gebv(gebv)

        # coerce the rest to thier correct shape
        method = self._coerce_1d_string_(method, "method", len(gebv))
        algorithm = self._coerce_1d_string_(algorithm, "algorithm", len(gebv))
        seed = self._coerce_1d_array(seed, "seed", len(gebv))
        cycle = self._coerce_1d_array(cycle, "cycle", len(gebv))
        score = self._coerce_2d_array(score, "score", len(gebv))
        if geno is not None:
            geno = self._coerce_geno(geno, len(gebv))

        self._selection_method.append(method)
        self._selection_algorithm.append(algorithm)
        self._selection_seed.append(seed)
        self._selection_cycle.append(cycle)
        self._selection_score.append(score)
        self._selection_gebv.append(gebv)
        if geno is not None:
            self._selection_genotype.append(geno)

    def concatenate_selection(self):
        if len(self._selection_method) > 1:
            self._selection_method = [numpy.concatenate(self._selection_method, axis=0)]
        if len(self._selection_algorithm) > 1:
            self._selection_algorithm = [numpy.concatenate(self._selection_algorithm, axis=0)]
        if len(self._selection_seed) > 1:
            self._selection_seed = [numpy.concatenate(self._selection_seed, axis=0)]
        if len(self._selection_cycle) > 1:
            self._selection_cycle = [numpy.concatenate(self._selection_cycle, axis=0)]
        if len(self._selection_score) > 1:
            self._selection_score = [numpy.concatenate(self._selection_score, axis=0)]
        if len(self._selection_gebv) > 1:
            self._selection_gebv = [numpy.concatenate(self._selection_gebv, axis=0)]
        if len(self._selection_genotype) > 1:
            self._selection_genotype = [numpy.concatenate(self._selection_genotype, axis=1)]

    def selection_history_to_dict(self, zfill = 6):
        self.concatenate_selection()

        # create selection stats dict
        selection_dict = {
            "method" : self._selection_method[0],
            "algorithm" : self._selection_algorithm[0],
            "seed" : self._selection_seed[0],
            "cycle" : self._selection_cycle[0]
        }

        # add scores column by column
        selection_dict.update(dict(zip(
            # column names
            ["score"+str(e+1) for e in range(self._selection_score[0].shape[1])],
            # columns
            self._selection_score[0].T
        )))

        # add GEBVs column by column
        selection_dict.update(dict(zip(
            # column names
            ["gebv"+str(e+1) for e in range(self._selection_gebv[0].shape[1])],
            # iterate through columns
            self._selection_gebv[0].T
        )))

        # add genotypes if we have them
        if len(self._selection_genotype) > 0:
            shape = self._selection_genotype[0].shape

            # add first phase keys to selection_dict
            selection_dict.update(dict(zip(
                # column names
                ["p1_"+str(e).zfill(zfill) for e in range(shape[1])],
                # iterate through columns
                self._selection_gebv[0][0].T
            )))

            # add second phase keys to selection_dict
            selection_dict.update(dict(zip(
                # column names
                ["p2_"+str(e).zfill(zfill) for e in range(shape[1])],
                # iterate through columns
                self._selection_gebv[0][1].T
            )))

        return selection_dict

    def selection_history_to_df(self, zfill = 6):
        # make dictionary
        selection_dict = self.selection_history_to_dict(zfill)

        # make pandas.DataFrame from dict
        selection_df = pandas.DataFrame(selection_dict)

        return selection_df

    def selection_history_to_csv(self, fname, zfill = 6, **kwargs):
        selection_df = self.selection_history_to_df(zfill)
        selection_df.to_csv(fname, **kwargs)

    ########################################################
    # general methods
    def reset(self):
        """
        Reset all internal history lists and arrays.
        """
        self.population_history_reset()
        self.selection_history_reset()

    def concatenate(self):
        """
        Concatenate internal history arrays for future exporting.
        """
        self.concatenate_population()
        self.concatenate_selection()

    def history_to_df(self, zfill = 6):
        population_df = self.population_history_to_df(zfill)
        selection_df = self.selection_history_to_df(zfill)

        return population_df, selection_df

    def history_to_csv(population_fname, selection_fname, zfill = 6, **kwargs):
        population_df, selection_df = self.history_to_df(zfill = zfill)
        population_df.to_csv(population_fname, **kwargs)
        selection_df.to_csv(selection_fname, **kwargs)

    def copy(self):
        """
        Copy the object.
        """
        cp = self.__class__(
            population = self._population,
            cross = self._cross,
            method = self._method
        )
        return cp

    def from_self(self, population = None, cross = None, method = None, **kwargs):
        deriv = self.__class__(
            population = self._population if population is None else population,
            cross = self._cross if cross is None else cross,
            method = self._method if method is None else method,
            **kwargs
        )
        return deriv

    def optimize(self, algorithm, **kwargs):
        """
        Optimize the breeding objective function provided an algorithm.

        Parameters
        ----------
        algorithm : Algorithm
            An optimizing algorithm.
        **kwargs
            Additional key word arguments for either the opimizer or the
            objective function.

        Returns
        -------
        sel : numpy.ndarray
            An array of selection decisions identified by the algorithm.

        Remarks
        -------
        Standard optimizing protocols. Override if necessary.
        """
        # ALWAYS begin with reseting the optimizing algorithm to remove history
        algorithm.reset()

        # we pass objcoeff onto optimizer. This will handle multiobjective.
        algorithm.optimize(
            self.objfn,
            **kwargs
        )

        # get global best selection configuration
        sel = algorithm.gbest_config()

        return sel

    def simulate(self, algorithm, bcycle = 0, objfn_kwargs = None, algo_kwargs = None, savegeno = False, seed = None):
        """
        Run breeding simulations.
        """
        # set seed if needed
        pybropt.util.cond_seed_rng(seed)

        # duplicate population and cross pointers
        population = self._population
        cross = self._cross

        # get initial population score
        pop_score = self.objfn(None, **objfn_kwargs)

        # get initial population GEBVs
        pop_gebv = population.gebv(None)

        # record history
        self.add_population_history(
            method = self._method,
            algorithm = algorithm.name,
            seed = seed,
            cycle = 0,
            score = pop_score,
            gebv = pop_gebv,
            geno = population.geno if savegeno else None
        )

        # simulate the breeding cycles
        for cycle in range(1, bcycle+1):
            # create a new Breeding object from self
            breeding = self.from_self(
                population = population,
                cross = cross
            )

            # make selections
            sel = breeding.optimize(
                algorithm = algorithm,
                **algo_kwargs,
                **objfn_kwargs
            )

            print("Selected:", sel)
            # get score of the selection
            sel_score = algorithm.gbest_score()

            # get selection GEBVs
            sel_gebv = self._population.gebv(sel)

            # add history selection
            self.add_selection_history(
                method = self._method,
                algorithm = algorithm.name,
                seed = seed,
                cycle = cycle,
                score = sel_score,
                gebv = sel_gebv,
                geno = population.geno[:,sel,:] if savegeno else None
            )

            # mate and generate new population
            population = cross.mate(sel)

            # get score of the entire population
            pop_score = self.objfn(None, **objfn_kwargs)

            # get entire population GEBVs
            pop_gebv = self._population.gebv(None)

            # record history
            self.add_population_history(
                method = self._method,
                algorithm = algorithm.name,
                seed = seed,
                cycle = cycle,
                score = pop_score,
                gebv = pop_gebv,
                geno = population.geno if savegeno else None
            )

            # make new cross object
            cross = cross.from_self(
                population = population,
            )
