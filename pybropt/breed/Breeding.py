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

    def method_population():
        doc = "The method_population property."
        def fget(self):
            return self._method_population
        def fset(self, value):
            self._method_population = value
        def fdel(self):
            del self._method_population
        return locals()
    method_population = property(**method_population())

    def method_selection():
        doc = "The method_selection property."
        def fget(self):
            return self._method_selection
        def fset(self, value):
            self._method_selection = value
        def fdel(self):
            del self._method_selection
        return locals()
    method_selection = property(**method_selection())

    def algorithm_population():
        doc = "The algorithm_population property."
        def fget(self):
            return self._algorithm_population
        def fset(self, value):
            self._algorithm_population = value
        def fdel(self):
            del self._algorithm_population
        return locals()
    algorithm_population = property(**algorithm_population())

    def algorithm_selection():
        doc = "The algorithm_selection property."
        def fget(self):
            return self._algorithm_selection
        def fset(self, value):
            self._algorithm_selection = value
        def fdel(self):
            del self._algorithm_selection
        return locals()
    algorithm_selection = property(**algorithm_selection())

    def seed_population():
        doc = "The seed_population property."
        def fget(self):
            return self._seed_population
        def fset(self, value):
            self._seed_population = value
        def fdel(self):
            del self._seed_population
        return locals()
    seed_population = property(**seed_population())

    def seed_selection():
        doc = "The seed_selection property."
        def fget(self):
            return self._seed_selection
        def fset(self, value):
            self._seed_selection = value
        def fdel(self):
            del self._seed_selection
        return locals()
    seed_selection = property(**seed_selection())

    def cycle_population():
        doc = "The cycle_population property."
        def fget(self):
            return self._cycle_population
        def fset(self, value):
            self._cycle_population = value
        def fdel(self):
            del self._cycle_population
        return locals()
    cycle_population = property(**cycle_population())

    def cycle_selection():
        doc = "The cycle_selection property."
        def fget(self):
            return self._cycle_selection
        def fset(self, value):
            self._cycle_selection = value
        def fdel(self):
            del self._cycle_selection
        return locals()
    cycle_selection = property(**cycle_selection())

    def score_population():
        doc = "The score_population property."
        def fget(self):
            return self._score_population
        def fset(self, value):
            self._score_population = value
        def fdel(self):
            del self._score_population
        return locals()
    score_population = property(**score_population())

    def score_selection():
        doc = "The score_selection property."
        def fget(self):
            return self._score_selection
        def fset(self, value):
            self._score_selection = value
        def fdel(self):
            del self._score_selection
        return locals()
    score_selection = property(**score_selection())

    def genotype_population():
        doc = "The genotype_population property."
        def fget(self):
            return self._genotype_population
        def fset(self, value):
            self._genotype_population = value
        def fdel(self):
            del self._genotype_population
        return locals()
    genotype_population = property(**genotype_population())

    def genotype_selection():
        doc = "The genotype_selection property."
        def fget(self):
            return self._genotype_selection
        def fset(self, value):
            self._genotype_selection = value
        def fdel(self):
            del self._genotype_selection
        return locals()
    genotype_selection = property(**genotype_selection())

    def gebv_population():
        doc = "The gebv_population property."
        def fget(self):
            return self._gebv_population
        def fset(self, value):
            self._gebv_population = value
        def fdel(self):
            del self._gebv_population
        return locals()
    gebv_population = property(**gebv_population())

    def gebv_selection():
        doc = "The gebv_selection property."
        def fget(self):
            return self._gebv_selection
        def fset(self, value):
            self._gebv_selection = value
        def fdel(self):
            del self._gebv_selection
        return locals()
    gebv_selection = property(**gebv_selection())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    def history_add_population(self, method, algorithm, seed, cycle, score,
        gebv, geno = None):
        # copy arrays
        method = numpy.string_(method)          # force byte string conversion
        algorithm = numpy.string_(algorithm)    # force byte string conversion
        seed = numpy.array(seed)
        cycle = numpy.array(cycle)
        score = numpy.array(score)
        gebv = numpy.array(gebv)
        if geno is not None:
            geno = numpy.array(geno)

        # adjustments for size and shape.
        if method.size == 1:
            method = numpy.repeat(method, len(gebv))
        if algorithm.size == 1:
            algorithm = numpy.repeat(algorithm, len(gebv))
        if seed.size == 1:
            seed = numpy.repeat(seed, len(gebv))
        if cycle.size == 1:
            cycle = numpy.repeat(cycle, len(gebv))
        if (geno is not None) and (geno.ndim == 2):
            geno = geno[:,None] # add axis to make 3d

        self._method_population.append(method)
        self._algorithm_population.append(algorithm)
        self._seed_population.append(seed)
        self._cycle_population.append(cycle)
        self._score_population.append(score)
        self._gebv_population.append(gebv)
        if geno is not None:
            self._genotype_population.append(geno)

    def history_add_selection(self, method, algorithm, seed, cycle, score,
        gebv, geno = None):
        # copy arrays
        method = numpy.string_(method)          # force byte string conversion
        algorithm = numpy.string_(algorithm)    # force byte string conversion
        seed = numpy.array(seed)
        cycle = numpy.array(cycle)
        score = numpy.array(score)
        gebv = numpy.array(gebv)
        if geno is not None:
            geno = numpy.array(geno)

        # adjustments for size and shape.
        if method.size == 1:
            method = numpy.repeat(method, len(gebv))
        if algorithm.size == 1:
            algorithm = numpy.repeat(algorithm, len(gebv))
        if seed.size == 1:
            seed = numpy.repeat(seed, len(gebv))
        if cycle.size == 1:
            cycle = numpy.repeat(cycle, len(gebv))
        if (geno is not None) and (geno.ndim == 2):
            geno = geno[:,None] # add axis to make 3d

        self._method_selection.append(method)
        self._algorithm_selection.append(algorithm)
        self._seed_selection.append(seed)
        self._cycle_selection.append(cycle)
        self._score_selection.append(score)
        self._gebv_selection.append(gebv)
        if geno is not None:
            self._genotype_selection.append(geno)

    def reset(self):
        """
        Reset all internal history lists and arrays.
        """
        self.reset_population()
        self.reset_selection()

    def reset_population(self):
        self._method_population = []
        self._algorithm_population = []
        self._seed_population = []
        self._cycle_population = []
        self._score_population = []
        self._genotype_population = []
        self._gebv_population = []

    def reset_selection(self):
        self._method_selection = []
        self._algorithm_selection = []
        self._seed_selection = []
        self._cycle_selection = []
        self._score_selection = []
        self._genotype_selection = []
        self._gebv_selection = []

    def concatenate(self):
        """
        Concatenate internal history arrays for future exporting.
        """
        if len(self._method_selection) > 1:
            self._method_selection = [numpy.concatenate(self._method_selection, axis=0)]
        if len(self._method_population) > 1:
            self._method_population = [numpy.concatenate(self._method_population, axis=0)]
        if len(self._algorithm_selection) > 1:
            self._algorithm_selection = [numpy.concatenate(self._algorithm_selection, axis=0)]
        if len(self._algorithm_population) > 1:
            self._algorithm_population = [numpy.concatenate(self._algorithm_population, axis=0)]
        if len(self._seed_selection) > 1:
            self._seed_selection = [numpy.concatenate(self._seed_selection, axis=0)]
        if len(self._seed_population) > 1:
            self._seed_population = [numpy.concatenate(self._seed_population, axis=0)]
        if len(self._cycle_selection) > 1:
            self._cycle_selection = [numpy.concatenate(self._cycle_selection, axis=0)]
        if len(self._cycle_population) > 1:
            self._cycle_population = [numpy.concatenate(self._cycle_population, axis=0)]
        if len(self._score_selection) > 1:
            self._score_selection = [numpy.concatenate(self._score_selection, axis=0)]
        if len(self._score_population) > 1:
            self._score_population = [numpy.concatenate(self._score_population, axis=0)]
        if len(self._genotype_selection) > 1:
            self._genotype_selection = [numpy.concatenate(self._genotype_selection, axis=1)]
        if len(self._genotype_population) > 1:
            self._genotype_population = [numpy.concatenate(self._genotype_population, axis=1)]

    def history_to_df(self, zfill = 6):
        # concatenate all histories
        self.concatenate()

        ####################################################
        # create population stats dict
        population_dict = {
            "method" : self._method_population[0],
            "algorithm" : self._algorithm_population[0],
            "seed" : self._seed_population[0],
            "cycle" : self._cycle_population[0]
        }

        # add scores to population_dict
        shape = self._score_population[0].shape
        if len(shape) == 1:
            # add a single key to population_dict
            population_dict.update({"score" : self._score_population[0]})
        elif len(shape) == 2:
            # add keys to population_dict
            population_dict.update(dict(zip(
                ["score" + e for e in range(shape[1])], # column names
                self._score_population[0].T             # iterate through columns
            )))

        # add GEBVs to population_dict
        shape = self._gebv_population[0].shape
        if len(shape) == 1:
            # add a single key to population_dict
            population_dict.update({"gebv" : self._gebv_population[0]})
        elif len(shape) == 2:
            # add keys to population_dict
            population_dict.update(dict(zip(
                ["gebv" + e for e in range(shape[1])],  # column names
                self._gebv_population[0].T              # iterate through columns
            )))

        # add genotypes to population_dict
        if len(self._genotype_population) > 0:
            shape = self._genotype_population[0].shape

            # add keys to population_dict
            population_dict.update(dict(zip(
                ["p1_"+str(e).zfill(zfill) for e in range(shape[1])],   # column names
                self._gebv_population[0][0].T                           # iterate through columns
            )))

            # add keys to population_dict
            population_dict.update(dict(zip(
                ["p2_"+str(e).zfill(zfill) for e in range(shape[1])],   # column names
                self._gebv_population[0][1].T                           # iterate through columns
            )))
        ####################################################


        ####################################################
        # create selection stats dict
        selection_dict = {
            "method" : self._method_selection[0],
            "algorithm" : self._algorithm_selection[0],
            "seed" : self._seed_selection[0],
            "cycle" : self._cycle_selection[0]
        }

        # add scores to selection_dict
        shape = self._score_selection[0].shape
        if len(shape) == 1:
            # add a single key to selection_dict
            selection_dict.update({"score" : self._score_selection[0]})
        elif len(shape) == 2:
            # add keys to selection_dict
            selection_dict.update(dict(zip(
                ["score" + e for e in range(shape[1])], # column names
                self._score_selection[0].T              # iterate through columns
            )))

        # add GEBVs to selection_dict
        shape = self._gebv_selection[0].shape
        if len(shape) == 1:
            # add a single key to selection_dict
            selection_dict.update({"gebv" : self._gebv_selection[0]})
        elif len(shape) == 2:
            # add keys to selection_dict
            selection_dict.update(dict(zip(
                ["gebv" + e for e in range(shape[1])],  # column names
                self._gebv_selection[0].T               # iterate through columns
            )))

        # add genotypes to selection_dict
        if len(self._genotype_selection) > 0:
            shape = self._genotype_selection[0].shape

            # add keys to population_dict
            selection_dict.update(dict(zip(
                ["p1_"+str(e).zfill(zfill) for e in range(shape[1])],   # column names
                self._gebv_selection[0][0].T                            # iterate through columns
            )))

            # add keys to population_dict
            selection_dict.update(dict(zip(
                ["p2_"+str(e).zfill(zfill) for e in range(shape[1])],   # column names
                self._gebv_selection[0][1].T                            # iterate through columns
            )))
        ####################################################

        # create DataFrame from dictionaries
        population_df = pandas.DataFrame(population_dict)
        selection_df = pandas.DataFrame(selection_dict)

        return population_df, selection_df

    def history_to_csv(population_fname, selection_fname, zfill = 6, *args, **kwargs):
        population_df, selection_df = self.history_to_df(zfill = zfill)

        population_df.to_csv(population_fname, *args, **kwargs)
        selection_df.to_csv(selection_fname, *args, **kwargs)

        return

    def objfn(self, sel, objcoeff = None, minimizing = True):
        """
        Breeding method objective function. Implement this in derived classes.

        Parameters
        ----------
        sel : numpy.ndarray
            Vector of selection data.
        *args : tuple
            Additional arguments for this function.
        **kwargs : dict
            Additional arguments for this function.

        Returns
        -------
        error : NotImplementedError
            This function raises an error on execution.
        """
        raise NotImplementedError("The method 'objfn' is abstract.")

    def objfn_vec(self, sel, objcoeff = None, minimizing = True):
        """
        Breeding method objective function. Implement this in derived classes.

        Parameters
        ----------
        sel : numpy.ndarray
            A 2D matrix of selection data.
        *args : tuple
            Additional arguments for this function.
        **kwargs : dict
            Additional arguments for this function.

        Returns
        -------
        error : NotImplementedError
            This function raises an error on execution.
        """
        raise NotImplementedError("The method 'objfn_vec' is abstract.")

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
        self.history_add_population(
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
            self.history_add_selection(
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
            self.history_add_population(
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
