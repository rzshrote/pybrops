class OPV(GenomicSelection):
    """docstring for OPV."""

    def __init__(self):
        super(OPV, self).__init__()

    def objfn(self, sel, geno = None, coeff = None, hcoeff = None):
        """
        Score a population of individuals based on Optimal Population Value (OPV)
        (Goiffon et al., 2017). Scoring for OPV is defined as the Genomic Estimated
        Breeding Value (GEBV) of the best inbred progeny produced from a breeding
        population allowing for an infinite number of intercrossings and
        generations. Haplotypes can be specified in OPV and markers within the
        haplotypes are assumed to always segregate with the haplotype. In simpler
        terms, OPV collects all of the best haplotype blocks and calculates the GEBV
        for those haplotypes for a given set of individuals.

        OPV selects the 'q' individuals to maximize the maximum GEBV possible within
        a population consisting of the 'q' selected individuals.

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
        hcoeff : numpy.ndarray
            A haplotype effect coefficients matrix of shape (t, m, n, p).
            Must be provided unless 'geno' and 'coeff' have been provided.
            Where:
                't' is the number of traits.
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.

        Returns
        -------
        opv : numpy.ndarray
            An OPV score matrix of shape (t,).
            Where:
                't' is the number of traits.
        """
        # if hcoeff is None, calculate it
        if hcoeff is None:
            hcoeff = OPV.hcoeff(geno, coeff)

        # if sel is None, slice all individuals
        if sel is None:
            sel = slice(None)

        # get the number of chromosome phases
        M = hcoeff.shape[1]

        # OPV calculation explanation
        # Step 1) find max value across all phases and individuals (axes=(1,2))
        # Step 2) sum across all loci (axis=1)
        # Step 3) multiply by number of phases (M)
        opv = M * numpy.amax(hcoeff[:,:,sel,:], (1,2)).sum(1)

        return opv

    def objfn_vec(self, sel, geno = None, coeff = None, hcoeff = None):
        """
        Score a population of individuals based on Optimal Population Value (OPV)
        (Goiffon et al., 2017). Scoring for OPV is defined as the Genomic Estimated
        Breeding Value (GEBV) of the best inbred progeny produced from a breeding
        population allowing for an infinite number of intercrossings and
        generations. Haplotypes can be specified in OPV and markers within the
        haplotypes are assumed to always segregate with the haplotype. In simpler
        terms, OPV collects all of the best haplotype blocks and calculates the GEBV
        for those haplotypes for a given set of individuals.

        OPV selects the 'q' individuals to maximize the maximum GEBV possible within
        a population consisting of the 'q' selected individuals.

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
        hcoeff : numpy.ndarray
            A haplotype effect coefficients matrix of shape (t, m, n, p).
            Must be provided unless 'geno' and 'coeff' have been provided.
            Where:
                't' is the number of traits.
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.

        Returns
        -------
        opv : numpy.ndarray
            An OPV score matrix of shape (j,t).
            Where:
                'j' is the number of selection configurations.
                't' is the number of traits.
        """
        # if hcoeff is None, calculate it
        if hcoeff is None:
            hcoeff = OPV.hcoeff(geno, coeff)

        # if sel is None, slice all individuals
        if sel is None:
            sel = numpy.arange(hcoeff.shape[2])[None,:]

        # get the number of chromosome phases
        M = hcoeff.shape[1]

        # OPV calculation explanation
        # Step 1) make selections
        # Step 2) find max value across all phases and individuals (axes=(1,3))
        # Step 3) sum across all loci (axis=2)
        # Step 4) transpose to get (j,t) shape
        # Step 4) multiply by number of phases (M)
        opv = M * numpy.amax(hcoeff[:,:,sel,:], (1,3)).sum(2).T

        return opv


    # TODO: add options to lump by specified haplotype groups
    @staticmethod
    def hcoeff(geno, coeff):
        """
        Calculate a haplotype coefficient matrix.

        Parameters
        ----------
        geno : numpy.ndarray
            A int8 binary genotype matrix of shape (m, n, p).
            Where:
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.
            Remarks:
                Shape of the matrix is most critical. Underlying matrix
                operations will support other numeric data types.
        coeff : numpy.ndarray
            A trait prediction coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.

        Returns
        -------
        hcoeff : numpy.ndarray
            A haplotype effect coefficients matrix of shape (t, m, n, p).
            Where:
                't' is the number of traits.
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.
        """
        # reshape matrices to shapes (1,m,n,p) and (t,1,1,p), respectively.
        hcoeff = geno[None,:] * coeff.T[:,None,None,:]

        # return product
        return hcoeff
