class MOGS(GenomicSelection):
    """docstring for MOGS."""

    def __init__(self):
        super(MOGS, self).__init__()

    def objfn(self, sel, geno, coeff = None, wcoeff = None, tfreq = None):
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
        mogs : numpy.ndarray
            A MOGS score matrix of shape (2,t)
        """
        # if no selection, select all
        if sel is None:
            sel = slice(None)

        # if wcoeff is None, calculate it
        if wcoeff is None:
            wcoeff = MOGS.wcoeff(coeff)

        # if wcoeff is None, calculate it
        if tfreq is None:
            tfreq = MOGS.tfreq(coeff)

        # generate a view of the geno matrix that only contains 'sel' rows.
        sgeno = geno[:,sel,:] # (m,n,p)[1] -> (m,k,p)

        # calculate number of phases
        phases = numpy.float64(sgeno.shape[0] * sgeno.shape[1])

        # calculate population frequencies; add axis for correct broadcast
        pfreq = (sgeno.sum((0,1)) / phases)[:,None] # (m,k,p)[0,1] -> (p,1)

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
        dist = numpy.absolute(tfreq - pfreq) # (p,t)-(p,1) -> (p,t)

        # compute f_PAU(x)
        pau = (wcoeff * allele_unavail).sum(0) # (p,t)[0] -> (t,)

        # compute f_PAFD(x)
        pafd = (wcoeff * dist).sum(0) # (p,t)[0] -> (t,)

        # stack to make MOGS matrix
        mogs = numpy.stack([pau, pafd]) # (2,t)

        return mogs

    def objfn_vec(self, sel, geno, coeff = None, wcoeff = None, tfreq = None):
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
        mogs : numpy.ndarray
            A MOGS score matrix of shape (2,j,t)
        """
        # if no selection, select all
        if sel is None:
            sel = numpy.arange(geno.shape[1])[None,:]

        # if wcoeff is None, calculate it
        if wcoeff is None:
            wcoeff = MOGS.wcoeff(coeff)

        # if wcoeff is None, calculate it
        if tfreq is None:
            tfreq = MOGS.tfreq(coeff)

        # generate a view of the geno matrix that only contains 'sel' rows.
        sgeno = geno[:,sel,:] # (m,n,p)[1] -> (m,j,k,p)

        # calculate number of phases
        phases = numpy.float64(sgeno.shape[0] * sgeno.shape[1])

        # calculate population frequencies; add axis for correct broadcast
        pfreq = (sgeno.sum((0,2)) / phases)[:,None] # (m,j,k,p)[0,2] -> (j,p,1)

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
        ) # (j,p,t)

        # calculate distance between target and population
        dist = numpy.absolute(tfreq - pfreq) # (p,t)-(j,p,1) -> (j,p,t)

        # compute f_PAU(x)
        pau = (wcoeff * allele_unavail).sum(1) # (j,p,t)[1] -> (j,t)

        # compute f_PAFD(x)
        pafd = (wcoeff * dist).sum(0) # (j,p,t)[1] -> (j,t)

        # stack to make MOGS matrix
        mogs = numpy.stack([pau, pafd]) # (2,j,t)

        return mogs


    @staticmethod
    def wcoeff(coeff):
        wcoeff = numpy.absolute(coeff)

        return wcoeff

    @staticmethod
    def tfreq(coeff):
        tfreq = (coeff >= 0).astype('float64')

        return tfreq