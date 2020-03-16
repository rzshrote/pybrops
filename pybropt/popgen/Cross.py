class Cross:
    """docstring for Cross."""
    ############################################################################
    ############################# Class Constants ##############################
    ############################################################################

    KEY_TO_VARAFN = {

    }

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################

    def __init__(self, population, varAfn = None, sparse = False,
        s = 0, t = 0, mem = None):
        """
        population : Population
            Population from which to determine
        """
        # check data types
        check_is_Population(population, "population")

        # set private variables
        self._population = population
        self._s = s
        self._t = t
        self._sparse = sparse
        if mem is None:
            self._mem = len(self._population.genetic_map)

    ############################################################################
    ############################ Object Properties #############################
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

    def variance():
        doc = "The variance property."
        def fget(self):
            return self._variance
        def fset(self, value):
            self._variance = value
        def fdel(self):
            del self._variance
        return locals()
    variance = property(**variance())

    def sparse():
        doc = "The sparse property."
        def fget(self):
            return self._sparse
        def fset(self, value):
            self._sparse = value
        def fdel(self):
            del self._sparse
        return locals()
    sparse = property(**sparse())

    def s():
        doc = "The s property."
        def fget(self):
            return self._s
        def fset(self, value):
            self._s = value
        def fdel(self):
            del self._s
        return locals()
    s = property(**s())

    def t():
        doc = "The t property."
        def fget(self):
            return self._t
        def fset(self, value):
            self._t = value
        def fdel(self):
            del self._t
        return locals()
    t = property(**t())

    def mem():
        doc = "The mem property."
        def fget(self):
            return self._mem
        def fset(self, value):
            self._mem = value
        def fdel(self):
            del self._mem
        return locals()
    mem = property(**mem())

    def varA_2wayDH_mat():
        doc = "The varA_2wayDH_mat property."
        def fget(self):
            if not hasattr(self, '_varA_2wayDH_mat'):
                self.calc_varA_2wayDH_mat()
            return self._varA_2wayDH_mat
        def fset(self, value):
            self._varA_2wayDH_mat = value
        def fdel(self):
            if hasattr(self, '_varA_2wayDH_mat'):
                del self._varA_2wayDH_mat
        return locals()
    varA_2wayDH_mat = property(**varA_2wayDH_mat())

    def varA_3wayDH_mat():
        doc = "The varA_3wayDH_mat property."
        def fget(self):
            if not hasattr(self, '_varA_3wayDH_mat'):
                self.calc_varA_3wayDH_mat()
            return self._varA_3wayDH_mat
        def fset(self, value):
            self._varA_3wayDH_mat = value
        def fdel(self):
            if hasattr(self, '_varA_3wayDH_mat'):
                del self._varA_3wayDH_mat
        return locals()
    varA_3wayDH_mat = property(**varA_3wayDH_mat())

    def varA_4wayDH_mat():
        doc = "The varA_4wayDH_mat property."
        def fget(self):
            if not hasattr(self, '_varA_4wayDH_mat'):
                self.calc_varA_4wayDH_mat()
            return self._varA_4wayDH_mat
        def fset(self, value):
            self._varA_4wayDH_mat = value
        def fdel(self):
            if hasattr(self, '_varA_4wayDH_mat'):
                del self._varA_4wayDH_mat
        return locals()
    varA_4wayDH_mat = property(**varA_4wayDH_mat())

    def varA_2wayDH_sparse_mat():
        doc = "The varA_2wayDH_sparse_mat property."
        def fget(self):
            if not hasattr(self, 'varA_')
            return self._varA_2wayDH_sparse_mat
        def fset(self, value):
            self._varA_2wayDH_sparse_mat = value
        def fdel(self):
            del self._varA_2wayDH_sparse_mat
        return locals()
    varA_2wayDH_sparse_mat = property(**varA_2wayDH_sparse_mat())

    def varA_3wayDH_sparse_mat():
        doc = "The varA_3wayDH_sparse_mat property."
        def fget(self):
            return self._varA_3wayDH_sparse_mat
        def fset(self, value):
            self._varA_3wayDH_sparse_mat = value
        def fdel(self):
            del self._varA_3wayDH_sparse_mat
        return locals()
    varA_3wayDH_sparse_mat = property(**varA_3wayDH_sparse_mat())

    def varA_4wayDH_sparse_mat():
        doc = "The varA_4wayDH_sparse_mat property."
        def fget(self):
            return self._varA_4wayDH_sparse_mat
        def fset(self, value):
            self._varA_4wayDH_sparse_mat = value
        def fdel(self):
            del self._varA_4wayDH_sparse_mat
        return locals()
    varA_4wayDH_sparse_mat = property(**varA_4wayDH_sparse_mat())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    @classmethod
    def varA(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def reset(self):
        """
        Remove all calculated matrices.
        """
        if hasattr(self, "varA_2wayDH_mat"):
            del self._varA_2wayDH_mat
        if hasattr(self, "varA_3wayDH_mat"):
            del self._varA_3wayDH_mat
        if hasattr(self, "varA_4wayDH_mat"):
            del self._varA_4wayDH_mat
        if hasattr(self, "varA_2wayDH_sparse_mat"):
            del self._varA_2wayDH_mat
        if hasattr(self, "varA_3wayDH_sparse_mat"):
            del self._varA_3wayDH_mat
        if hasattr(self, "varA_4wayDH_sparse_mat"):
            del self._varA_4wayDH_mat

    @classmethod
    def calc_varA_2wayDH_mat(self):
        """
        Calculate a symmetrical matrix of progeny variance for each pairwise
        2-way cross. Calculations are derived from Osthushenrich et al. (2017).

        Assumes inbred individuals.
        """
        # create aliases for several variables
        geno = self._population.geno
        coeff = self._population.genomic_model.coeff
        gmap = self._population.genetic_map
        mem = self._mem

        # get the number of individuals
        ntaxa = geno.shape[1]

        # allocate a square matrix for each pairwise variance
        self._varA_2wayDH_mat = numpy.zeros((ntaxa, ntaxa), dtype='float64')

        # for each linkage group
        for lst, lsp in zip(gmap.chr_grp_stix, gmap.chr_grp_spix):
            # for each computational chunk
            for rst,rsp in zip(range(lst,lsp,mem),srange(lst+mem,lsp,mem)):
                for cst,csp in zip(range(lst,lsp,mem),srange(lst+mem,lsp,mem)):
                    # get recombination probability matrix for chunk
                    r = gmap.recomb_prob(rst, rsp, cst, csp)

                    # calculate a D1 matrix; this is specific to mating scheme
                    D1 = D_1(r, self._s, self._t)

                    # get marker coefficients for rows and columns
                    rcoeff = coeff[rst:rsp]
                    ccoeff = coeff[cst:csp]

                    # for each mate pair (excluding selfs)
                    for female in range(1,n_indiv): # varA row index
                        for male in range(0,female): # varA col index
                            # calculate genotype differences for row, col {-1,0,1}
                            rdgeno = geno[0,female,rst:rsp] - geno[0,male,rst:rsp]
                            cdgeno = geno[0,female,cst:csp] - geno[0,male,cst:csp]

                            # calculate effect differences
                            reffect = rdgeno * rcoeff
                            ceffect = cdgeno * ccoeff

                            # compute the dot products to get a partial variance
                            varA_part = reffect.dot(D).dot(ceffect)

                            # add this partial variance to the lower triangle
                            self._varA_2wayDH_mat[female,male] += varA_part

        # since varA matrix is symmetrical, copy lower triangle to the upper
        for female in range(1, n_indiv):
            for male in range(0, female):
                self._varA_2wayDH_mat[male,female] = self._varA_2wayDH_mat[female,male]

        return

    @classmethod
    def calc_varA_3wayDH_mat(self):
        """
        Calculate a symmetrical matrix of progeny variance for each pairwise
        3-way cross. Calculations are derived from Osthushenrich et al. (2017).

        Assumes inbred individuals.
        """
        # create aliases for several variables
        geno = self._population.geno
        coeff = self._population.genomic_model.coeff
        gmap = self._population.genetic_map
        mem = self._mem

        # get the number of individuals
        ntaxa = geno.shape[1]

        # allocate a square matrix for each pairwise variance
        self._varA_3wayDH_mat = numpy.zeros(
            (ntaxa, ntaxa, ntaxa),
            dtype='float64'
        )

        # for each linkage group
        for lst, lsp in zip(gmap.chr_grp_stix, gmap.chr_grp_spix):
            # for each computational chunk
            for rst,rsp in zip(range(lst,lsp,mem),srange(lst+mem,lsp,mem)):
                for cst,csp in zip(range(lst,lsp,mem),srange(lst+mem,lsp,mem)):
                    # get recombination probability matrix for chunk
                    r = gmap.recomb_prob(rst, rsp, cst, csp)

                    # calculate a D1 matrix; this is specific to mating scheme
                    D1 = D_1(r, self._s, self._t)

                    # calculate a D2 matrix; this is specific to mating scheme
                    D2 = D_2(r, self._s, self._t)

                    # get marker coefficients for rows and columns
                    rcoeff = coeff[rst:rsp]
                    ccoeff = coeff[cst:csp]

                    # for each 3-way cross (excluding selfs)
                    # subscript codes:
                    #   1 = recurrent parent
                    #   2 = female
                    #   3 = male
                    for recurr in range(0,ntaxa):           # varA slice index
                        for female in range(0,ntaxa):       # varA row index
                            # calculate genotype differences for row, col
                            rdgeno21 = geno[0,female,rst:rsp] - geno[0,recurr,rst:rsp]
                            cdgeno21 = geno[0,female,cst:csp] - geno[0,recurr,cst:csp]

                            # calculate effect differences
                            reffect21 = rdgeno13 * rcoeff
                            ceffect21 = cdgeno13 * ccoeff

                            # calculate varA part for female-recurrent
                            varA_part21 = reffect21.dot(D1).dot(ceffect21)

                            # only do lower triangle since symmetrical within each slice
                            for male in range(0,female):    # varA col index
                                # calculate genotype differences for row, col
                                rdgeno23 = geno[0,female,rst:rsp] - geno[0,male,rst:rsp]
                                cdgeno23 = geno[0,female,cst:csp] - geno[0,male,cst:csp]

                                rdgeno31 = geno[0,male,rst:rsp] - geno[0,recurr,rst:rsp]
                                cdgeno31 = geno[0,male,cst:csp] - geno[0,recurr,cst:csp]

                                # calculate effect differences
                                reffect23 = rdgeno12 * rcoeff
                                ceffect23 = cdgeno12 * ccoeff

                                reffect31 = rdgeno23 * rcoeff
                                ceffect31 = cdgeno23 * ccoeff

                                # calculate varA parts for crosses with male
                                varA_part23 = reffect23.dot(D2).dot(ceffect23)
                                varA_part31 = reffect31.dot(D1).dot(ceffect31)

                                # calculate varA part for this matrix chunk
                                varA_part = (2.0 * (varA_part21 + varA_part31)) + varA_part23

                                # add this partial variance to the lower triangle
                                self._varA_3wayDH_mat[recurr,female,male] += varA_part

        # divide entire matrix by 4 to get variance per the equation
        self._varA_3wayDH_mat /= 4.0

        # each matrix is symmetrical within a slice because exchanging female
        # and male orders is mathematically equivalent.
        # copy lower triangle to the upper since varA matrix is symmetrical within each slice
        for female in range(1, ntaxa):
            for male in range(0, female):
                self._varA_3wayDH_mat[:,male,female] = self._varA_3wayDH_mat[:,female,male]

        return

    @classmethod
    def calc_varA_4wayDH_mat(self):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def calc_varA_2wayDH_sparse_mat(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def calc_varA_3wayDH_sparse_mat(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def calc_varA_4wayDH_sparse_mat(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_2wayDH(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_3wayDH(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_4wayDH(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_2wayDH_sparse(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_3wayDH_sparse(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_4wayDH_sparse(self, sel):
        raise NotImplementedError("Method not implemented.")



    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################

    @staticmethod
    def r_k(r, k):
        """
        Calculate the expected observed recombination rate between two loci at
        filial generation 'k' as specified by Lehermeier (2017).

        Parameters
        ----------
        r : numpy.ndarray
            Recombination probability matrix.
        k : int, inf
            Selfing filial generation number to derive gametes from.
            Example:
                k = 1    ->  Derive gametes from F1
                k = 2    ->  Derive gametes from F2
                ...
                k = inf  ->  Derive gametes from SSD

        Returns
        -------
        r_k : numpy.ndarray
            Expected observed recombination rate between two loci at filial
            generation 'k'.
        """
        # multiply matrix by two
        two_r = 2.0 * r

        # calculate first component of r_k
        r_k = two_r / (1.0 + two_r)

        # if k != inf, we do not have SSD and second term is needed
        if k != numpy.inf:
            r_k *= (1.0 - ( (0.5**k) * ((1.0 - two_r)**k) ) )

        # return result
        return r_k

    def D_1(r, s, t):
        """
        Calculate a D_1 matrix. This matrix represents a linkage disequilibrium
        prototype for one recombination event prior to selfing or random
        intermating.

        Parameters
        ----------
        r : numpy.ndarray
            Recombination probability matrix.
        s : int, inf
            Selfing generation number to derive gametes from.
            Example:
                s = 0    ->  Derive gametes from F1
                s = 1    ->  Derive gametes from F2
                s = 2    ->  Derive gametes from F3
                ...
                k = inf  ->  Derive gametes from SSD
            Remark:
                's' takes priority over 't'. If s > 0, 't' calculations are not
                performed.
        t : int, inf
            Random intermating generation number to derive gametes from.
            Example:
                t = 0    ->  Derive gametes from F1
                t = 1    ->  Derive gametes from (F1 x F1)
                t = 2    ->  Derive gametes from (F1 x F1) x (F1 x F1)
                ...
                t = inf  ->  Derive gametes from unlimited random intermating
            Remark:
                's' takes priority over 't'. If s > 0, 't' calculations are not
                performed.

        Returns
        -------
        D1 : numpy.ndarray
            A D_1 matrix.
        """
        if s == 0 and t == 0:
            return (1 - 2*r)
        elif s > 0:
            rk = Cross.r_k(r, s+1)
            return (1.0 - 2.0*rk)
        elif t > 0:
            return (1.0 - 2.0*r) * ((1.0 - r)**t)
        else:
            raise ValueError("s and t must be > 0.")

    def D_2(r, s, t):
        """
        Calculate a D_2 matrix. This matrix represents a linkage disequilibrium
        prototype for two recombination events prior to selfing or random
        intermating.

        Parameters
        ----------
        r : numpy.ndarray
            Recombination probability matrix.
        s : int, inf
            Selfing generation number to derive gametes from.
            Example:
                s = 0    ->  Derive gametes from F1
                s = 1    ->  Derive gametes from F2
                s = 2    ->  Derive gametes from F3
                ...
                k = inf  ->  Derive gametes from SSD
            Remark:
                's' takes priority over 't'. If s > 0, 't' calculations are not
                performed.
        t : int, inf
            Random intermating generation number to derive gametes from.
            Example:
                t = 0    ->  Derive gametes from F1
                t = 1    ->  Derive gametes from (F1 x F1)
                t = 2    ->  Derive gametes from (F1 x F1) x (F1 x F1)
                ...
                t = inf  ->  Derive gametes from unlimited random intermating
            Remark:
                's' takes priority over 't'. If s > 0, 't' calculations are not
                performed.

        Returns
        -------
        D2 : numpy.ndarray
            A D_2 matrix.
        """
        if s == 0 and t == 0:
            return (1.0 - 2.0*r)**2
        elif s > 0:
            rk = Cross.r_k(r, s+1)
            four_r = 4.0*r
            return ( 1.0 - four_r + (four_r*rk) )
        elif t > 0:
            return ((1.0 - 2.0*r)**2) * ((1.0 - r)**t)
        else:
            raise ValueError("s and t must be > 0.")
