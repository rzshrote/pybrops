class PAU(GenomicSelection):
    """docstring for PAU."""

    def __init__(self):
        super(PAU, self).__init__()

    def objfn(self, sel, geno, coeff = None, wcoeff = None, tfreq = None):
        """
        Population Allele Unvailability (PAU) objective function.
            The goal is to minimize this function. Lower is better.
            This is a bare bones function. Minimal error checking.

        Given the provided genotype matrix 'geno' and row selections from it 'sel',
        calculate the selection allele freq. From the selection allele frequencies
        and the target allele frequencies, determine if the target frequencies
        cannot be attained after unlimited generations and selection rounds.
        Multiply this vector by a weight coefficients vector 'wcoeff'.

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
        wcoeff : numpy.ndarray
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
        pau : numpy.ndarray
            A PAU score matrix of shape (t,)
        """
        # if weight coefficient is not provided, calculate it.
        if wcoeff is None:
            wcoeff = PAU.wcoeff(coeff)

        # if 'sel' is None, set to all individuals
        if sel is None:
            sel = slice(None)

        # generate a view of the geno matrix that only contains 'sel' rows.
        sgeno = geno[:,sel,:]

        # calculate num chromosomes in the selection (m * n) as double
        phases = numpy.float64(sgeno.shape[0] * sgeno.shape[1])

        # calculate population frequencies; add axis for correct broadcast
        pfreq = (sgeno.sum((0,1)) / phases)[:,None]

        # calculate some inequalities for use multiple times
        pfreq_lteq_0 = (pfreq <= 0.0)   # is population freq <= 0.0
        pfreq_gteq_1 = (pfreq >= 1.0)   # is population freq >= 1.0

        # calculate allele unavailability
        allele_unavail = numpy.where(
            tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
            pfreq_lteq_0,           # then set True if sel has allele freq == 0
            numpy.where(            # else
                tfreq > 0.0,        # if 0.0 < target freq < 1.0
                numpy.logical_or(   # then set TRUE if pop freq is outside (0.0,1.0)
                    pfreq_lteq_0,
                    pfreq_gteq_1
                ),
                pfreq_gteq_1        # else set TRUE if pop freq is >= 1.0
            )
        )

        # calculate PAU score
        pau = (wcoeff * allele_unavail).sum(0)

        # return score as the specified output data type
        return pau

    def objfn_vec(self, sel, geno, coeff = None, wcoeff = None, tfreq = None):
        """
        Population Allele Unvailability (PAU) objective function.
            The goal is to minimize this function. Lower is better.
            This is a bare bones function. Minimal error checking.

        Given the provided genotype matrix 'geno' and row selections from it 'sel',
        calculate the selection allele freq. From the selection allele frequencies
        and the target allele frequencies, determine if the target frequencies
        cannot be attained after unlimited generations and selection rounds.
        Multiply this vector by a weight coefficients vector 'wcoeff'.

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
        wcoeff : numpy.ndarray
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
        pau : numpy.ndarray
            A PAU score matrix of shape (t,)
        """
        # if 'sel' is None, set to all individuals
        if sel is None:
            sel = numpy.arange(geno.shape[1])[None,:]

        # if weight coefficient is not provided, calculate it.
        if wcoeff is None:
            wcoeff = PAU.wcoeff(coeff)

        # if tfreq is not provided, calculate it.
        if tfreq is None:
            tfreq = PAU.tfreq(coeff)

        # generate a view of the geno matrix that only contains 'sel' rows.
        sgeno = geno[:,sel,:]

        # calculate num chromosomes in the selection (m * n) as double
        phases = numpy.float64(sgeno.shape[0] * sgeno.shape[2])

        # calculate population frequencies; add axis for correct broadcast
        pfreq = (sgeno.sum((0,2)) / phases)[:,None]

        # calculate some inequalities for use multiple times
        pfreq_lteq_0 = (pfreq <= 0.0)   # is population freq <= 0.0
        pfreq_gteq_1 = (pfreq >= 1.0)   # is population freq >= 1.0

        # calculate allele unavailability
        allele_unavail = numpy.where(
            tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
            pfreq_lteq_0,           # then set True if sel has allele freq == 0
            numpy.where(            # else
                tfreq > 0.0,        # if 0.0 < target freq < 1.0
                numpy.logical_or(   # then set TRUE if pop freq is outside (0.0,1.0)
                    pfreq_lteq_0,
                    pfreq_gteq_1
                ),
                pfreq_gteq_1        # else set TRUE if pop freq is >= 1.0
            )
        )

        # calculate PAU score
        pau = (wcoeff * allele_unavail).sum(1)

        # return score as the specified output data type
        return pau

    @staticmethod
    def wcoeff(coeff):
        wcoeff = numpy.absolute(coeff)

        return wcoeff

    @staticmethod
    def tfreq(coeff):
        tfreq = (coeff >= 0).astype('float64')

        return tfreq
