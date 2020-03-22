class PAFD(GenomicSelection):
    """docstring for PAFD."""

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################
    @classmethod
    def __init__(self, population, cross, wcoeff = None, tfreq = None):
        super(PAFD, self).__init__(population, cross)

        # check that we have marker coefficients
        check_is_ParametricGenomicModel(self._population.genomic_model)

        # set wcoeff if needed
        if wcoeff is None:
            wcoeff = self.calc_wcoeff()
        self._wcoeff = wcoeff

        # set tfreq if needed
        if tfreq is None:
            tfreq = self.calc_tfreq()
        self._tfreq = tfreq

    ############################################################################
    ################################ Properties ################################
    ############################################################################

    def wcoeff():
        doc = "The wcoeff property."
        def fget(self):
            return self._wcoeff
        def fset(self, value):
            self._wcoeff = value
        def fdel(self):
            del self._wcoeff
        return locals()
    wcoeff = property(**wcoeff())

    def tfreq():
        doc = "The tfreq property."
        def fget(self):
            return self._tfreq
        def fset(self, value):
            self._tfreq = value
        def fdel(self):
            del self._tfreq
        return locals()
    tfreq = property(**tfreq())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    @classmethod
    def calc_wcoeff(self):
        """
        Calculate weight coefficients.
        """
        # calculate wcoeff
        wcoeff = PAFD.wcoeff_mat(self._population.coeff)

        return wcoeff

    @classmethod
    def calc_tfreq(self):
        """
        Calculate target frequencies
        """
        # calculate tfreq
        tfreq = PAFD.tfreq_mat(self._population.coeff)

        return tfreq

    @classmethod
    def objfn(self, sel, objcoeff = None, negate = True):
        # calculate PAFD values
        pafd = PAFD.objfn_mat(
            sel,
            self._population.geno,
            wcoeff = self._wcoeff,
            tfreq = self._tfreq
        )

        # negate PAFD scores if necessary.
        if negate:
            pafd = -pafd

        # if we have objective weights, take dot product for weight sum method
        if objcoeff is not None:
            pafd = pafd.dot(objcoeff)

        return pafd

    @classmethod
    def objfn_vec(self, sel, objcoeff = None, negate = True):
        # calculate PAFD values
        pafd = PAFD.objfn_vec_mat(
            sel,
            self._population.geno,
            wcoeff = self._wcoeff,
            tfreq = self._tfreq
        )

        # negate OPV scores if necessary
        if negate:
            pafd = -pafd

        # take the dot product if necessary
        if objcoeff is not None:
            pafd = pafd.dot(objcoeff)

        return pafd

    @classmethod
    def optimize(self, objcoeff = None, negate = True, algorithm = None,
        gbestix = None, *args, **kwargs):
        # we pass objcoeff onto optimizer. This will handle multiobjective.
        algorithm.optimize(
            self.objfn,
            *args,
            **kwargs,
            objcoeff = objcoeff,
            negate = negate
        )

        # get global best
        gbest = algorithm.gbest()

        # get selection indices or whole tuple
        sel = gbest[gbestix] if gbestix is not None else gbest

        return sel

    @classmethod
    def simulate(self):
        raise NotImplementedError

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def objfn_mat(sel, geno, coeff = None, wcoeff = None, tfreq = None):
        """
        Population Allele Frequency Distance (PAFD) objective function.

        The goal is to minimize this function. Lower is better.

        This is a bare bones function. Minimal error checking.

        Parameters
        ==========
        sel : numpy.ndarray, None
            A selection indices matrix of shape (k,)
            Where:
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, use all individuals.
        geno : numpy.ndarray, None
            A int8 binary genotype matrix of shape (m, n, p).
            Where:
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.
            Remarks:
                Shape of the matrix is most critical. Underlying matrix
                operations will support other numeric data types.
        coeff : numpy.ndarray, None
            A trait prediction coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        wcoeff : numpy.ndarray, None
            A marker weight coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            Remarks: Values in 'wcoeff' have an assumption:
                All values must be non-negative.
        tfreq : None, floating, numpy.ndarray
            A target allele frequency matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            Example:
                tfreq = numpy.array([0.2, 0.6, 0.7])

        Returns
        =======
        pafd : numpy.ndarray
            A PAFD score matrix of shape (t,)
        """
        # if no selection, select all
        if sel is None:
            sel = slice(None)

        # if wcoeff is None, calculate it
        if wcoeff is None:
            wcoeff = PAFD.wcoeff_mat(coeff)

        # if tfreq is None, calculate it
        if tfreq is None:
            tfreq = PAFD.tfreq_mat(coeff)

        # generate a view of the geno matrix that only contains 'sel' rows.
        sgeno = geno[:,sel,:]

        # calculate the number of chromosome phases present (depth * rows)
        phases = numpy.float64(sgeno.shape[0] * sgeno.shape[1])

        # calculate population frequencies; add axis for correct broadcast
        pfreq = (sgeno.sum((0,1)) / phases)[:,None] # (p,1)

        # calculate distance between target and population
        dist = numpy.absolute(tfreq - pfreq) # (p,t)-(p,1) -> (p,t)

        # calculate PAFD score
        pafd = (wcoeff * dist).sum(0) # (p,t)[0] -> (t,)

        # return score as the specified output data type
        return pafd

    @staticmethod
    def objfn_vec_mat(sel, geno, coeff = None, wcoeff = None, tfreq = None):
        """
        Population Allele Frequency Distance (PAFD) objective function.

        The goal is to minimize this function. Lower is better.

        This is a bare bones function. Minimal error checking.

        Parameters
        ==========
        sel : numpy.ndarray, None
            A selection indices matrix of shape (j,k)
            Where:
                'j' is the number of selection configurations.
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, use all individuals.
        geno : numpy.ndarray, None
            A int8 binary genotype matrix of shape (m, n, p).
            Where:
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.
            Remarks:
                Shape of the matrix is most critical. Underlying matrix
                operations will support other numeric data types.
        coeff : numpy.ndarray, None
            A trait prediction coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        wcoeff : numpy.ndarray, None
            A marker weight coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            Remarks: Values in 'wcoeff' have an assumption:
                All values must be non-negative.
        tfreq : None, floating, numpy.ndarray
            A target allele frequency matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            Example:
                tfreq = numpy.array([0.2, 0.6, 0.7])

        Returns
        =======
        pafd : numpy.ndarray
            A PAFD score matrix of shape (j,t)
        """
        # if wcoeff is None, calculate it
        if wcoeff is None:
            wcoeff = PAFD.wcoeff_mat(coeff)

        # if tfreq is not provided, calculate it.
        if tfreq is None:
            tfreq = PAFD.tfreq_mat(coeff)

        # generate a view of the geno matrix that only contains 'sel' rows.
        sgeno = geno[:,sel,:] # (m,n,p)[1] -> (m,j,k,p)

        # calculate the number of chromosome phases present (depth * rows)
        phases = numpy.float64(sgeno.shape[0] * sgeno.shape[2])

        # calculate population frequencies; add axis for correct broadcast
        pfreq = (sgeno.sum((0,2)) / phases)[:,None] # (m,j,k,p)[0,2] -> (j,p,1)

        # calculate distance between target and population
        dist = numpy.absolute(tfreq - pfreq) # (p,t)-(j,p,1) -> (j,p,t)

        # calculate PAFD score
        pafd = (wcoeff * dist).sum(1) # (j,p,t)[1] -> (j,t)

        # return score as the specified output data type
        return pafd

    @staticmethod
    def wcoeff_mat(coeff):
        """
        Calculate weight coefficients matrix from a coefficients matrix by
        taking the absolute value of the elements in the matrix.
        """
        wcoeff = numpy.absolute(coeff)

        return wcoeff

    @staticmethod
    def tfreq_mat(coeff):
        """
        Calculate target frequencies.
        """
        tfreq = (coeff >= 0).astype('float64')

        return tfreq
