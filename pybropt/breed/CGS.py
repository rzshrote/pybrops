class CGS(GenomicSelection):
    """docstring for CGS."""

    def __init__(self):
        super(CGS, self).__init__()

    def objfn(self, sel, geno, coeff):
        """
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
        Genomic Estimated Breeding Values (GEBV) for a population.

        CGS selects the 'q' individuals with the largest GEBVs.

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

        Returns
        -------
        cgs : numpy.ndarray
            Returns an array of floating point number representing GEBVs for each
            individual.
        """
        # if sel is None, slice all individuals
        if sel is None:
            sel = slice(None)
            #return geno.sum(0).dot(coeff)

        # CGS calculation explanation
        # Step 1: get sum of alleles for each individual: shape=(len(sel),L)
        # Step 2: take dot product to get CGS values
        cgs = geno[:,sel].sum(0).dot(coeff)

        return cgs

    def objfn_vec(self, sel, geno, coeff):
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
        geno : numpy.ndarray
            An array of allele states. The dtype of 'geno' should be 'uint8'. Array
            shape should be (depth, row, column) = (M, N, L) where 'M' represents
            number of chromosome phases, 'N' represents number of individuals, 'L'
            represents number of markers. Array format should be the 'C' format.
        coeff : numpy.ndarray
            An array of coefficients for allele effects. The dtype of 'coeff'
            should be either 'float32' or 'float64'. This array should be single
            dimensional (column,) = (N,) where 'N' represents number of markers.

        Returns
        -------
        cgs : numpy.ndarray
            Returns an array of floating point number representing GEBVs for each
            individual.
        """
        if sel is None:
            return geno.sum(0).dot(coeff)
        return geno[:,sel].sum(0).dot(coeff)

    def optimize(self, algorithm):
        pass
