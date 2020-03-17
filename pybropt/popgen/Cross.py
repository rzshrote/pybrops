import numpy

class Cross:
    """docstring for Cross."""
    ############################################################################
    ############################# Class Constants ##############################
    ############################################################################

    KEY_TO_VARAFN = {
        '2wayDH' :      ('varA_2wayDH', 'varA_2wayDH_vec'),
        '3wayDH' :      ('varA_3wayDH', 'varA_3wayDH_vec'),
        '4wayDH' :      ('varA_4wayDH', 'varA_4wayDH_vec'),
        'dihybridDH' :  ('varA_dihybridDH', 'varA_dihybridDH_vec'),
        None :          ('varA_default', 'varA_default')
    }

    KEY_TO_VARAFN_SPARSE = {
        '2wayDH' :      ('varA_2wayDH_sparse', 'varA_2wayDH_sparse_vec'),
        '3wayDH' :      ('varA_3wayDH_sparse', 'varA_3wayDH_sparse_vec'),
        '4wayDH' :      ('varA_4wayDH_sparse', 'varA_4wayDH_sparse_vec'),
        'dihybridDH' :  ('varA_dihybridDH_sparse', 'varA_dihybridDH_sparse_vec'),
        None :          ('varA_default', 'varA_default')
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
        self.set_varAfn(varAfn, sparse)
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

    def varA_dihybridDH_mat():
        doc = "The varA_dihybridDH_mat property."
        def fget(self):
            if not hasattr(self, "_varA_dihybridDH_mat"):
                self.calc_varA_dihybridDH_mat()
            return self._varA_dihybridDH_mat
        def fset(self, value):
            self._varA_dihybridDH_mat = value
        def fdel(self):
            del self._varA_dihybridDH_mat
        return locals()
    varA_dihybridDH_mat = property(**varA_dihybridDH_mat())

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

    def varA_dihybridDH_sparse_mat():
        doc = "The varA_dihybridDH_sparse_mat property."
        def fget(self):
            return self._varA_dihybridDH_sparse_mat
        def fset(self, value):
            self._varA_dihybridDH_sparse_mat = value
        def fdel(self):
            del self._varA_dihybridDH_sparse_mat
        return locals()
    varA_dihybridDH_sparse_mat = property(**varA_dihybridDH_sparse_mat())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    @classmethod
    def reset(self):
        """
        Remove all calculated matrices.
        """
        # remove matrices if they exist.
        if hasattr(self, "_varA_2wayDH_mat"):
            del self._varA_2wayDH_mat
        if hasattr(self, "_varA_3wayDH_mat"):
            del self._varA_3wayDH_mat
        if hasattr(self, "_varA_4wayDH_mat"):
            del self._varA_4wayDH_mat
        if hasattr(self, "_varA_dihybridDH_mat"):
            del self._varA_dihybridDH_mat
        if hasattr(self, "_varA_2wayDH_sparse_mat"):
            del self._varA_2wayDH_mat
        if hasattr(self, "_varA_3wayDH_sparse_mat"):
            del self._varA_3wayDH_mat
        if hasattr(self, "_varA_4wayDH_sparse_mat"):
            del self._varA_4wayDH_mat
        if hasattr(self, "_varA_dihybridDH_sparse_mat"):
            del self._varA_dihybridDH_sparse_mat

    ############################################################################
    ###################### Variance Related Class Methods ######################
    @classmethod
    def set_varAfn(self, varAfn = None, sparse = False):
        """
        Set the varA() and varA_vec() attributes. Calling these will result in
        the gathering of a default variance component.

        Parameters
        ----------
        varAfn : str, None
            A key to lookup in internal lookup tables.
            Valid keys: {'2wayDH', '3wayDH', '4wayDH', 'dihybridDH'}
        sparse : boolean, default = False
            Whether 'varAfn' is a sparse matrix calculation.
        """
        # make string variables to hold the name of attributes
        varA_str = None
        varA_vec_str = None

        # lookup varAfn in lookup table
        if sparse:
            varA_str, varA_vec_str = Cross.KEY_TO_VARAFN_SPARSE[varAfn]
        else:
            varA_str, varA_vec_str = Cross.KEY_TO_VARAFN[varAfn]

        # set new function attributes
        self.varA = getattr(self, varA_str)
        self.varA_vec = getattr(self, varA_vec_str)

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
                    D1_mat = Cross.D1(r, self._s, self._t)

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

        return self._varA_2wayDH_mat

    @classmethod
    def calc_varA_3wayDH_mat(self):
        """
        Calculate a symmetrical matrix of progeny variance for each pairwise
        3-way cross.

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
                    D1_mat = Cross.D1(r, self._s, self._t)

                    # calculate a D2 matrix; this is specific to mating scheme
                    D2_mat = Cross.D2(r, self._s, self._t)

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
                            reffect21 = rdgeno21 * rcoeff
                            ceffect21 = cdgeno21 * ccoeff

                            # calculate varA part for female-recurrent
                            varA_part21 = reffect21.dot(D1_mat).dot(ceffect21)

                            # only do lower triangle since symmetrical within each slice
                            for male in range(0,female):    # varA col index
                                # calculate genotype differences for row, col
                                rdgeno23 = geno[0,female,rst:rsp] - geno[0,male,rst:rsp]
                                cdgeno23 = geno[0,female,cst:csp] - geno[0,male,cst:csp]

                                rdgeno31 = geno[0,male,rst:rsp] - geno[0,recurr,rst:rsp]
                                cdgeno31 = geno[0,male,cst:csp] - geno[0,recurr,cst:csp]

                                # calculate effect differences
                                reffect23 = rdgeno23 * rcoeff
                                ceffect23 = cdgeno23 * ccoeff

                                reffect31 = rdgeno31 * rcoeff
                                ceffect31 = cdgeno31 * ccoeff

                                # calculate varA parts for crosses with male
                                varA_part23 = reffect23.dot(D2_mat).dot(ceffect23)
                                varA_part31 = reffect31.dot(D1_mat).dot(ceffect31)

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

        return self._varA_3wayDH_mat

    @classmethod
    def calc_varA_4wayDH_mat(self):
        """
        Calculate a symmetrical matrix of progeny variance for each pairwise
        4-way cross.

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
        self._varA_4wayDH_mat = numpy.zeros(
            (ntaxa, ntaxa, ntaxa, ntaxa),
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
                    D1_mat = Cross.D1(r, self._s, self._t)

                    # calculate a D2 matrix; this is specific to mating scheme
                    D2_mat = Cross.D2(r, self._s, self._t)

                    # get marker coefficients for rows and columns
                    rcoeff = coeff[rst:rsp]
                    ccoeff = coeff[cst:csp]

                    # for each 4-way cross (excluding selfs)
                    # subscript codes:
                    #   1 = female 2
                    #   2 = male 2
                    #   3 = female 1
                    #   4 = male 1
                    # TODO: make this more computationally efficient.
                    #       operations can be simplified if you mirror sections
                    #       of the 4d matrix between blocks.
                    #       how to do this is difficult to conceive
                    for female2 in range(0,ntaxa):              # varA block index (change me for efficiency?)
                        for male2 in range(0,ntaxa):            # varA slice index (change me for efficiency?)
                            # calculate genotype differences for row, col
                            rdgeno21 = geno[0,male2,rst:rsp] - geno[0,female2,rst:rsp]
                            cdgeno21 = geno[0,male2,cst:csp] - geno[0,female2,cst:csp]

                            # calculate effect differences
                            reffect21 = rdgeno21 * rcoeff
                            ceffect21 = cdgeno21 * ccoeff

                            # calculate varA part for male2-female2
                            varA_part21 = reffect21.dot(D2_mat).dot(ceffect21)

                            for female1 in range(0,ntaxa):      # varA row index
                                # calculate genotype differences for row, col
                                rdgeno31 = geno[0,female1,rst:rsp] - geno[0,female2,rst:rsp]
                                cdgeno31 = geno[0,female1,cst:csp] - geno[0,female2,cst:csp]
                                rdgeno32 = geno[0,female1,rst:rsp] - geno[0,male2,rst:rsp]
                                cdgeno32 = geno[0,female1,cst:csp] - geno[0,male2,cst:csp]

                                # calculate effect differences
                                reffect31 = rdgeno31 * rcoeff
                                ceffect31 = cdgeno31 * ccoeff
                                reffect32 = rdgeno32 * rcoeff
                                ceffect32 = cdgeno32 * ccoeff

                                # calculate varA part for female1
                                varA_part31 = reffect31.dot(D1_mat).dot(ceffect31)
                                varA_part32 = reffect32.dot(D1_mat).dot(ceffect32)

                                # only do lower triangle since symmetrical within each slice
                                for male1 in range(0,female1):  # varA col index
                                    # calculate genotype differences for row, col
                                    rdgeno41 = geno[0,male1,rst:rsp] - geno[0,female2,rst:rsp]
                                    cdgeno41 = geno[0,male1,cst:csp] - geno[0,female2,cst:csp]
                                    rdgeno42 = geno[0,male1,rst:rsp] - geno[0,male2,rst:rsp]
                                    cdgeno42 = geno[0,male1,cst:csp] - geno[0,male2,cst:csp]
                                    rdgeno43 = geno[0,male1,rst:rsp] - geno[0,female1,rst:rsp]
                                    cdgeno43 = geno[0,male1,cst:csp] - geno[0,female1,cst:csp]

                                    # calculate effect differences
                                    reffect41 = rdgeno41 * rcoeff
                                    ceffect41 = rdgeno41 * ccoeff
                                    reffect42 = rdgeno42 * rcoeff
                                    ceffect42 = rdgeno42 * ccoeff
                                    reffect43 = rdgeno43 * rcoeff
                                    ceffect43 = rdgeno43 * ccoeff

                                    # calculate varA parts for crosses with male
                                    varA_part41 = reffect41.dot(D1_mat).dot(ceffect41)
                                    varA_part42 = reffect42.dot(D1_mat).dot(ceffect42)
                                    varA_part43 = reffect43.dot(D2_mat).dot(ceffect43)

                                    # calculate varA part for this matrix chunk
                                    varA_part = varA_part21 + varA_part31 + varA_part32 + varA_part41 + varA_part42 + varA_part43

                                    # add this partial variance to the lower triangle
                                    self._varA_4wayDH_mat[female2,male2,female1,male1] += varA_part

        # divide entire matrix by 4 to get variance per the equation
        self._varA_4wayDH_mat /= 4.0

        # each matrix is symmetrical within a slice because exchanging female
        # and male orders is mathematically equivalent.
        # copy lower triangle to the upper since varA matrix is symmetrical within each slice
        for female1 in range(1, ntaxa):
            for male1 in range(0, female1):
                self._varA_4wayDH_mat[:,:,male1,female1] = self._varA_4wayDH_mat[:,:,female1,male1]

        return self._varA_4wayDH_mat

    @classmethod
    def calc_varA_dihybridDH_mat(self):
        """
        Calculate a symmetrical matrix of progeny variance for each pairwise
        dihybrid cross.

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
        self._varA_dihybridDH_mat = numpy.zeros(
            (ntaxa, ntaxa),
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
                    D1_mat = Cross.D1(r, self._s, self._t)

                    # calculate a D2 matrix; this is specific to mating scheme
                    D2_mat = Cross.D2(r, self._s, self._t)

                    # get marker coefficients for rows and columns
                    rcoeff = coeff[rst:rsp]
                    ccoeff = coeff[cst:csp]

                    # for each mate pair (including selfs)
                    # subscript codes:
                    #   1 = female phase 2
                    #   2 = female phase 1
                    #   3 = male phase 2
                    #   4 = male phase 1
                    for female in range(0,n_indiv):     # varA row index
                        for male in range(0,female):    # varA col index
                            # calculate genotype differences for row, col
                            rdgeno21 = geno[0,female,rst:rsp] - geno[1,female,rst:rsp]
                            cdgeno21 = geno[0,female,cst:csp] - geno[1,female,cst:csp]
                            rdgeno31 = geno[1,male,rst:rsp] - geno[1,female,rst:rsp]
                            cdgeno31 = geno[1,male,cst:csp] - geno[1,female,cst:csp]
                            rdgeno32 = geno[1,male,rst:rsp] - geno[0,female,rst:rsp]
                            cdgeno32 = geno[1,male,cst:csp] - geno[0,female,cst:csp]
                            rdgeno41 = geno[0,male,rst:rsp] - geno[1,female,rst:rsp]
                            cdgeno41 = geno[0,male,cst:csp] - geno[1,female,cst:csp]
                            rdgeno42 = geno[0,male,rst:rsp] - geno[0,female,rst:rsp]
                            cdgeno42 = geno[0,male,cst:csp] - geno[0,female,cst:csp]
                            rdgeno43 = geno[0,male,rst:rsp] - geno[1,male,rst:rsp]
                            cdgeno43 = geno[0,male,cst:csp] - geno[1,male,cst:csp]

                            # calculate effect differences
                            reffect21 = rdgeno21 * rcoeff
                            ceffect21 = cdgeno21 * ccoeff
                            reffect31 = rdgeno31 * rcoeff
                            ceffect31 = cdgeno31 * ccoeff
                            reffect32 = rdgeno32 * rcoeff
                            ceffect32 = cdgeno32 * ccoeff
                            reffect41 = rdgeno41 * rcoeff
                            ceffect41 = rdgeno41 * ccoeff
                            reffect42 = rdgeno42 * rcoeff
                            ceffect42 = rdgeno42 * ccoeff
                            reffect43 = rdgeno43 * rcoeff
                            ceffect43 = rdgeno43 * ccoeff

                            # calculate varA part for female2-male
                            varA_part21 = reffect21.dot(D2_mat).dot(ceffect21)
                            varA_part31 = reffect31.dot(D1_mat).dot(ceffect31)
                            varA_part32 = reffect32.dot(D1_mat).dot(ceffect32)
                            varA_part41 = reffect41.dot(D1_mat).dot(ceffect41)
                            varA_part42 = reffect42.dot(D1_mat).dot(ceffect42)
                            varA_part43 = reffect43.dot(D2_mat).dot(ceffect43)

                            # calculate varA part for this matrix chunk
                            varA_part = varA_part21 + varA_part31 + varA_part32 + varA_part41 + varA_part42 + varA_part43

                            # add this partial variance to the lower triangle
                            self._varA_dyhybridDH_mat[female,male] += varA_part

        # divide entire matrix by 4 to get variance per the equation
        self._varA_dyhybridDH_mat /= 4.0

        # each matrix is symmetrical within a slice because exchanging female
        # and male orders is mathematically equivalent.
        # copy lower triangle to the upper since varA matrix is symmetrical within each slice
        for female in range(1, ntaxa):
            for male in range(0, female):
                self._varA_dyhybridDH_mat[male,female] = self._varA_dyhybridDH_mat[female,male]

        return self._varA_4wayDH_mat

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
    def calc_varA_dihybridDH_sparse_mat(self, sel):
        raise NotImplementedError("Method not implemented.")

    ##############################################
    @classmethod
    def varA_default(self, sel):
        """
        Default method for varA method. Raises a RuntimeError.
        """
        raise RuntimeError("varA function not set.")
    ##############################################

    ##############################################
    ########### Non-sparse matrix ops ############
    @classmethod
    def varA_2wayDH(self, sel):
        """
        Retrieve additive variance components for a 2-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7

        Returns
        -------
        varA_val : numpy.ndarray
            A 1D array of variance components of shape (k/2,)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_2wayDH_mat"):
            self.calc_varA_2wayDH_mat()

        # get female, male indices
        female = sel[0::2]
        male = sel[1::2]

        # get variance terms
        varA_val = self._varA_2wayDH_mat[female,male]

        # return variance terms
        return varA_val

    @classmethod
    def varA_3wayDH(self, sel):
        """
        Retrieve additive variance components for a 3-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the recurrent parent.
                Second index is the female parent.
                Third index is the male parent.
            Example:
                [1,5,3,8,2,7]
                recurrent = 1,8
                female = 5,2
                male = 3,7

        Returns
        -------
        varA_val : numpy.ndarray
            A 1D array of variance components of shape (k/3,)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_3wayDH_mat"):
            self.calc_varA_3wayDH_mat()

        # get recurrent, female, male indices
        recurr = sel[0::3]
        female = sel[1::3]
        male = sel[2::3]

        # get variance terms
        varA_val = self._varA_3wayDH_mat[recurr,female,male]

        # return variance terms
        return varA_val

    @classmethod
    def varA_4wayDH(self, sel):
        """
        Retrieve additive variance components for a 4-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of shape (k,).
            Where:
                'k' is the number of selected individuals, divisible by 4.
            Indices are paired as follows:
                First index is the female parent 2.
                Second index is the male parent 2.
                Third index is the female parent 1.
                Fourth index is the male parent 1.
            Example:
                [1,5,3,8]
                female2 = 1
                male2 = 5
                female1 = 3
                male1 = 8

        Returns
        -------
        varA_val : numpy.ndarray
            A 1D array of variance components of shape (k/4,)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_4wayDH_mat"):
            self.calc_varA_4wayDH_mat()

        # get female2, male2, female1, male1 indices
        female2 = sel[0::4]
        male2 = sel[1::4]
        female = sel[2::4]
        male = sel[3::4]

        # get variance terms
        varA_val = self._varA_3wayDH_mat[female2,male2,female1,male1]

        # return variance terms
        return varA_val

    @classmethod
    def varA_dihybridDH(self, sel):
        """
        Retrieve additive variance components for a dihybrid DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7

        Returns
        -------
        varA_val : numpy.ndarray
            A 1D array of variance components of shape (k/2,)
        """
        if not hasattr(self, "_varA_dihybridDH_mat"):
            self.calc_varA_dihybridDH_mat()

        # get female, male indices
        female = sel[0::2]
        male = sel[1::2]

        # get variance terms
        varA_val = self._varA_dihybridDH_mat[female,male]

        # return variance terms
        return varA_val

    @classmethod
    def varA_2wayDH_vec(self, sel):
        """
        Retrieve additive variance components for a 2-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 2D array of indices of shape (j,k).
            Where:
                'j' is the number of mating configurations.
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [[1,5,3,8,2,7],
                 [2,7,3,5,4,6]]
                female = [1,3,2],[2,3,4]
                male = [5,8,7],[7,5,6]

        Returns
        -------
        varA_val : numpy.ndarray
            A 2D array of variance components of shape (j,k/2)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_2wayDH_mat"):
            self.calc_varA_2wayDH_mat()

        # get female, male indices
        female = sel[:,0::2]
        male = sel[:,1::2]

        # get variance terms
        varA_val = self._varA_2wayDH_mat[female,male]

        # return variance terms
        return varA_val

    @classmethod
    def varA_3wayDH_vec(self, sel):
        """
        Retrieve additive variance components for a 3-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 2D array of indices of shape (j,k).
            Where:
                'j' is the number of mating configurations.
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the recurrent parent.
                Second index is the female parent.
                Third index is the male parent.
            Example:
                [[1,5,3,8,2,7],
                 [2,7,3,5,4,6]]
                recurrent = [1,8],[2,5]
                female = [5,2],[7,4]
                male = [3,7],[3,6]

        Returns
        -------
        varA_val : numpy.ndarray
            A 2D array of variance components of shape (j,k/3)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_3wayDH_mat"):
            self.calc_varA_3wayDH_mat()

        # get recurrent, female, male indices
        recurr = sel[:,0::3]
        female = sel[:,1::3]
        male = sel[:,2::3]

        # get variance terms
        varA_val = self._varA_3wayDH_mat[recurr,female,male]

        # return variance terms
        return varA_val

    @classmethod
    def varA_4wayDH_vec(self, sel):
        """
        Retrieve additive variance components for a 4-way DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 2D array of indices of shape (j,k).
            Where:
                'j' is the number of mating configurations.
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the female parent 2.
                Second index is the male parent 2.
                Third index is the female parent 1.
                Fourth index is the male parent 1.
            Example:
                [[1,5,3,8]
                 [7,2,9,0]]
                female2 = [1],[7]
                male2 = [5],[2]
                female1 = [3],[9]
                male1 = [8],[0]

        Returns
        -------
        varA_val : numpy.ndarray
            A 2D array of variance components of shape (j,k/4)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_4wayDH_mat"):
            self.calc_varA_4wayDH_mat()

        # get female2, male2, female1, male1 indices
        female2 = sel[:,0::4]
        male2 = sel[:,1::4]
        female = sel[:,2::4]
        male = sel[:,3::4]

        # get variance terms
        varA_val = self._varA_3wayDH_mat[female2,male2,female1,male1]

        # return variance terms
        return varA_val

    @classmethod
    def varA_dihybridDH_vec(self, sel):
        """
        Retrieve additive variance components for a dihybrid DH cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 2D array of indices of shape (j,k).
            Where:
                'j' is the number of mating configurations.
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [[1,5,3,8,2,7],
                 [2,7,3,5,4,6]]
                female = [1,3,2],[2,3,4]
                male = [5,8,7],[7,5,6]

        Returns
        -------
        varA_val : numpy.ndarray
            A 2D array of variance components of shape (j,k/2)
        """
        # if the matrix has not been calculated, calculate it.
        if not hasattr(self, "_varA_dihybridDH_mat"):
            self.calc_varA_dihybridDH_mat()

        # get female, male indices
        female = sel[:,0::2]
        male = sel[:,1::2]

        # get variance terms
        varA_val = self._varA_dihybridDH_mat[female,male]

        # return variance terms
        return varA_val
    ##############################################

    ##############################################
    ############# Sparse matrix ops ##############
    @classmethod
    def varA_2wayDH_sparse(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_3wayDH_sparse(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_4wayDH_sparse(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_dihybridDH_sparse(self, self):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_2wayDH_sparse_vec(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_3wayDH_sparse_vec(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_4wayDH_sparse_vec(self, sel):
        raise NotImplementedError("Method not implemented.")

    @classmethod
    def varA_dihybridDH_sparse_vec(self, self):
        raise NotImplementedError("Method not implemented.")
    ##############################################

    ############################################################################


    ############################################################################
    ####################### Mating Related Class Methods #######################
    @classmethod
    def meiosis(self, sel, seed = None):
        """
        Simulate meiosis. Generate gametes from a genotype matrix.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of shape (k,).
            Where:
                'k' is the number of selected individuals.

        Returns
        -------
        gametes : numpy.ndarray
        """
        # if there is a RNG seed, seed the RNG
        if seed is not None:
            numpy.random.seed(seed)

        # aliases for variables
        geno = self._population.geno
        gmap = self._population.genetic_map

        # make an empty matrix to contain gametes
        gamete = numpy.empty(
            (len(sel), geno.shape[2]),
            dtype = 'int8'
        )

        # make a vector of probabilities of recombination between neighbor loci
        gmap_stix = self._population.genetic_map.chr_grp_stix
        gmap_spix = self._population.genetic_map.chr_grp_spix
        gmap_grp_len = self._population.genetic_map.chr_grp_len
        gmap_mapfn = self._population.genetic_map.mapfn
        gmap_pos = self._population.genetic_map.map_pos

        prob_grp_len = gmap_grp_len - 1
        prob_spix = prob_grp_len.cumsum()
        prob_stix = prob_spix - prob_grp_len

        prob = numpy.empty(prob_grp_len.sum(), dtype = 'float64')

        # calculate probabilities for neighbor recombination
        for gst,gsp,pst,psp in zip(gmap_stix, gmap_spix, prob_stix, prob_spix):
            prob[pst:psp] = gmap_mapfn(gmap_pos[gst+1:gsp] - gmap_pos[gst:gsp-1])

        # generate random numbers
        rnd = numpy.random.uniform(0, 1, (len(sel),len(prob)))

        # for each gamete from a parent
        for i,s in enumerate(sel):
            # for each linkage group
            for gst,gsp,pst,psp in zip(gmap_stix, gmap_spix, prob_stix, prob_spix):
                # determine where crossovers occur (add 1 for exclusive index)
                xo = numpy.flatnonzero(rnd[i,pst:psp] <= prob[pst:psp]) + 1

                phase = numpy.random.binomial(1, 0.5) # starting phase
                spt = 0 # start point offset
                for ept in xo:
                    gamete[i,gst+spt:gst+ept] = geno[phase,s,gst+spt:gst+ept]
                    phase = 1 - phase # alternate phase
                    spt = ept # move to the next copying segment
                gamete[i,gst+spt:gsp] = geno[phase,s,gst+spt:gsp]

        return gamete

    @classmethod
    def mate_2way_controlled(self, sel, n, seed = None):
        """
        Perform a controlled cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7

        n : numpy.ndarray, int
            An array of shape (k/2,) designating the number of progeny to
            generate per cross. If an integer is provided, the number of each
            cross is equivalent.

        Returns
        -------
        population : Population
            A new population of individuals from the designated crosses.
        """
        # if a seed is provided, seed the RNG
        if seed is not None:
            numpy.random.seed(seed)

        # double the length of 'n'
        n = numpy.repeat(n, 2)

        # calculate gamete sources
        gametesrc = numpy.tile(sel, n)

        # generate male and female gametes.
        gamete = self.meiosis(gametesrc)

        # stack our gametes to progeny
        progeny = numpy.stack((gamete[0::2,:], gamete[1::2,:]))

        # create a new population
        population = Population(
            progeny,
            self._population.genomic_model,
            self._population.genetic_map
        )

        # return our progeny
        return population

    @classmethod
    def mate_DH(self, sel, n, s = 0, t = 0, seed = None):
        """
        Perform a controlled cross.

        Parameters
        ----------
        sel : numpy.ndarray
            A 1D array of indices of shape (k,).
            Where:
                'k' is the number of selected individuals.
        n : numpy.ndarray, int
            An array of shape (k/2,) designating the number of progeny to
            generate per cross. If an integer is provided, the number of each
            cross is equivalent.

        Returns
        -------
        population : Population
            A new population of individuals from the designated crosses.
        """
        raise NotImplementedError()
        # if a seed is provided, seed the RNG
        if seed is not None:
            numpy.random.seed(seed)

        # double the length of 'n'
        n = numpy.repeat(n, 2)

        # calculate gamete sources
        gametesrc = numpy.tile(sel, n)

        # generate male and female gametes.
        gamete = self.meiosis(gametesrc)

        # stack our gametes to progeny
        progeny = numpy.stack((gamete[0::2,:], gamete[1::2,:]))

        # create a new population
        population = Population(
            progeny,
            self._population.genomic_model,
            self._population.genetic_map
        )

        # return our progeny
        return population

    @classmethod
    def mate_2way_wrand(self, sel, n, weight, seed = None):
        """
        """
        raise NotImplementedError()

    @classmethod
    def mate_2way_erand(self, sel, n, exact, seed = None):
        """
        """
        # if a seed is provided, seed the RNG
        if seed is not None:
            numpy.random.seed(seed)

        if exact[0::2].sum() != exact[1::2].sum():
            raise ValueError("Gamete numbers do not match.")

        # get source indices
        female = numpy.repeat(sel[0::2], exact[0::2])
        male = numpy.repeat(sel[1::2], exact[1::2])

        # shuffle source indices
        numpy.random.shuffle(female)
        numpy.random.shuffle(male)

        # allocate mate config array
        mate_config = numpy.empty(len(female) + len(male), dtype = sel.dtype)

        # assign mate configurations
        mate_config[0::2] = female
        mate_config[1::2] = male

        # make new population
        population = self.mate_2way_controlled(mate_config, n)

        return population

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################

    @staticmethod
    def rk(r, k):
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

    def D1(r, s, t):
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
        D_1 : numpy.ndarray
            A D_1 matrix.
        """
        if s == 0 and t == 0:
            return (1 - 2*r)
        elif s > 0:
            rk = Cross.rk(r, s+1)
            return (1.0 - 2.0*rk)
        elif t > 0:
            return (1.0 - 2.0*r) * ((1.0 - r)**t)
        else:
            raise ValueError("s and t must be > 0.")

    def D2(r, s, t):
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
        D_2 : numpy.ndarray
            A D_2 matrix.
        """
        if s == 0 and t == 0:
            return (1.0 - 2.0*r)**2
        elif s > 0:
            rk = Cross.rk(r, s+1)
            four_r = 4.0*r
            return ( 1.0 - four_r + (four_r*rk) )
        elif t > 0:
            return ((1.0 - 2.0*r)**2) * ((1.0 - r)**t)
        else:
            raise ValueError("s and t must be > 0.")
