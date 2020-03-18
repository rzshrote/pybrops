class Breeding:
    """docstring for Breeding."""

    ############################################################################
    ############################# Reserved methods #############################
    ############################################################################

    def __init__(self):
        """
        Constructor for Breeding class.
        """
        self.reset()

    ############################################################################
    ################################ Properties ################################
    ############################################################################

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

    @classmethod
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

    @classmethod
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


    @classmethod
    def reset(self):
        """
        Reset all internal history lists and arrays.
        """
        self._method_selection = []
        self._method_population = []
        self._algorithm_selection = []
        self._algorithm_population = []
        self._seed_selection = []
        self._seed_population = []
        self._cycle_selection = []
        self._cycle_population = []
        self._score_selection = []
        self._score_population = []
        self._genotype_selection = []
        self._genotype_population = []

    @classmethod
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

    @classmethod
    def history_to_df(self, zfill = 6):
        # concatenate all histories
        self.concatenate()

        population_dict = {
            "method" : self._method_population[0],
            "algorithm" : self._algorithm_population[0],
            "seed" : self._seed_population[0],
            "cycle" : self._cycle_population[0],
            "score" : self._score_population[0],
            "gebv" : self._gebv_population[0]
        }
        selection_dict = {
            "method" : self._sel_method[0],
            "algorithm" : self._sel_algorithm[0],
            "seed" : self._sel_seed[0],
            "cycle" : self._sel_cycle[0],
            "score" : self._sel_score[0],
            "gebv" : self._sel_gebv[0]
        }

        if len(self._genotype_population) > 0:
            population_p1head = [
                "p1"+str(i).zfill(zfill)
                for i in range(self._pop_phase1[0].shape[1])
            ]
            population_p2head = [
                "p2"+str(i).zfill(zfill)
                for i in range(self._pop_phase2[0].shape[1])
            ]
            for i,header in enumerate(pop_p1head):
                population_dict[header] = self._genotype_population[0][0,:,i]
            for i,header in enumerate(pop_p2head):
                population_dict[header] = self._genotype_population[0][1,:,i]

        if len(self._genotype_selection) > 0:
            selection_p1head = [
                "p1"+str(i).zfill(zfill)
                for i in range(self._sel_phase1[0].shape[1])
            ]
            selection_p2head = [
                "p2"+str(i).zfill(zfill)
                for i in range(self._sel_phase2[0].shape[1])
            ]
            for i,header in enumerate(sel_p1head):
                selection_dict[header] = self._genotype_selection[0][0,:,i]
            for i,header in enumerate(sel_p2head):
                selection_dict[header] = self._genotype_selection[0][1,:,i]

        population_df = pandas.DataFrame(population_dict)
        selection_df = pandas.DataFrame(selection_dict)

        return population_df, selection_df

    @classmethod
    def history_to_csv(population_fname, selection_fname, zfill = 6):
        population_df, selection_df = self.history_to_df(zfill = zfill)

        population_df.to_csv(population_fname)
        selection_df.to_csv(selection_fname)

        return

    def objfn(self, sel, *args, **kwargs):
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
        raise NotImplementedError("This method is not implemented.")

    def objfn_vec(self, sel, *args, **kwargs):
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
        raise NotImplementedError("This method is not implemented.")

    @classmethod
    def optimize(self, algorithm):
        """
        """
        raise NotImplementedError("This method is not implemented.")

    @classmethod
    def simulate(self, algorithm):
        """
        """
        raise NotImplementedError("This method is not implemented.")
