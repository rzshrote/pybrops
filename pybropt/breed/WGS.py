class WGS(GenomicSelection):
    """docstring for CGS."""

    def __init__(self):
        super(CGS, self).__init__()

    def objfn(self, sel, geno, coeff = None, weight = None, wcoeff = None):
        """
        Score a population of individuals based on Weighted Genomic Selection (WGS)
        (Goddard, 2009; Jannink, 2010). Scoring for WGS is defined as the adjusted
        Genomic Estimated Breeding Values (GEBV) for an individual. Marker
        coefficients are adjusted according to the frequency of the most beneficial
        allele.

        WGS selects the 'q' individuals with the largest weighted GEBVs.

        Parameters
        ----------
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
        weight : numpy.ndarray, None
            A allele effect weight coefficients matrix of shape (p,t)
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        wcoeff : numpy.ndarray, None
            A allele weight coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            This is the Hadamard product of 'coeff' and 1/sqrt('weight')

        Returns
        -------
        wgs : numpy.ndarray
            A weighted GEBV matrix of shape (k,t)
            Returns an array of floating point number representing GEBVs for each
            individual.
        """
        # if sel is None, slice all individuals
        if sel is None:
            sel = slice(None)

        # if no wcoeff matrix, calculate it
        if wcoeff is None:
            wcoeff = WGS.wcoeff(coeff, weight)

        # CGS calculation explanation
        # Step 1: get sum of alleles for each individual: shape=(len(sel),L)
        # Step 2: take dot product to get CGS values
        wgs = geno[:,sel].sum(0).dot(wcoeff)

        return wgs

    def objfn_vec(self, sel, geno, coeff = None, weight = None, wcoeff = None):
        """
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
        Genomic Estimated Breeding Values (GEBV) for a population.

        CGS selects the 'q' individuals with the largest GEBVs.

        Parameters
        ----------
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
        weight : numpy.ndarray, None
            A allele effect weight coefficients matrix of shape (p,t)
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        wcoeff : numpy.ndarray, None
            A allele weight coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            This is the Hadamard product of 'coeff' and 'wcoeff'

        Returns
        -------
        wgs : numpy.ndarray
            Returns an array of floating point number representing GEBVs for each
            individual.
        """
        # if sel is None, slice all individuals
        if sel is None:
            sel = slice(None)

        # if no wcoeff matrix, calculate it
        if wcoeff is None:
            wcoeff = WGS.wcoeff(coeff, weight)

        # CGS calculation explanation
        # Step 1: get sum of alleles for each individual: shape=(len(sel),L)
        # Step 2: take dot product to get CGS values
        wgs = geno[:,sel].sum(0).dot(wcoeff)

        return wgs

    @staticmethod
    def wcoeff(coeff, weight):
        """
        Parameters
        ----------
        coeff : numpy.ndarray, None
            A trait prediction coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        weight : numpy.ndarray, None
            A allele effect weight coefficients matrix of shape (p,t)
            Where:
                'p' is the number of markers.
                't' is the number of traits.

        Returns
        -------
        wcoeff : numpy.ndarray
            A allele weight coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            This is the Hadamard product of 'coeff' and 'wcoeff'
        """
        # divide coeff by sqrt(weight)
        wcoeff = coeff / numpy.sqrt(weight)

        return weight

    def optimize(self, algorithm):
        pass
