class PAFD(GenomicSelection):
    """docstring for PAFD."""

    def __init__(self):
        super(PAFD, self).__init__()

    @classmethod
    def objfn(self, sel, geno, coeff = None, wcoeff = None, tfreq = None):
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
            wcoeff = PAFD.wcoeff(coeff)

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


    @classmethod
    def objfn_vec(self, sel, geno, coeff = None, wcoeff = None, tfreq = None):
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
            wcoeff = PAFD.wcoeff(coeff)

        # if tfreq is not provided, calculate it.
        if tfreq is None:
            tfreq = PAFD.tfreq(coeff)

        # generate a view of the geno matrix that only contains 'sel' rows.
        sgeno = geno[:,sel,:] # (m,n,p)[1] -> (m,j,k,p)

        # calculate the number of chromosome phases present (depth * rows)
        phases = numpy.float64(sgeno.shape[0] * sgeno.shape[1])

        # calculate population frequencies; add axis for correct broadcast
        pfreq = (sgeno.sum((0,2)) / phases)[:,None] # (m,j,k,p)[0,2] -> (j,p,1)

        # calculate distance between target and population
        dist = numpy.absolute(tfreq - pfreq) # (p,t)-(j,p,1) -> (j,p,t)

        # calculate PAFD score
        pafd = (wcoeff * dist).sum(1) # (j,p,t)[1] -> (j,t)

        # return score as the specified output data type
        return pafd

    @staticmethod
    def wcoeff(coeff):
        wcoeff = numpy.absolute(coeff)

        return wcoeff

    @staticmethod
    def tfreq(coeff):
        tfreq = (coeff >= 0).astype('float64')

        return tfreq
