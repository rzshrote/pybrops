from . import TwoWayDHCross

class GenerationalTwoWayDHCross(TwoWayDHCross):
    """docstring for GenerationalTwoWayDHCross."""

    def __init__(self, gmult, **kwargs):
        """
        Parameters
        ----------
        gmult : int
            Generation multiplier to apply to group labels
        """
        super(GenerationalTwoWayDHCross, self).__init__()
        self.gmult = gmult

    def mate(self, t_cur, t_max, pgvmat, sel, ncross, nprogeny, rng, s = 0, **kwargs):
        """
        Mate individuals according to a 2-way mate selection scheme.

        Parameters
        ----------
        pgvmat : PhasedGenotypeVariantMatrix
            A GenotypeVariantMatrix containing candidate breeding individuals.
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                Even indices are female.
                Odd indices are male.
            Example:
                [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7
        ncross : numpy.ndarray
            Number of cross patterns to perform.
        nprogeny : numpy.ndarray
            Number of doubled haploid progeny to generate per cross.
        s : int, default = 0
            Number of selfing generations post-cross before double haploids are
            generated.

        Returns
        -------
        progeny : PhasedGenotypeVariantMatrix
            A PhasedGenotypeVariantMatrix of progeny.
        """
        progeny = super(GenerationalTwoWayDHCross, self).mate(
            t_cur = t_cur,
            t_max = t_max,
            pgvmat = pgvmat,
            sel = sel,
            ncross = ncross,
            nprogeny = nprogeny,
            rng = rng,
            s = s
        )

        # add taxa_grp to progeny
        taxa_grp = numpy.repeat(numpy.arange(numpy.sum(ncross, dtype='int64')), nprogeny)
        taxa_grp += t_cur * self.gmult
        progeny.taxa_grp = taxa_grp

        return progeny
