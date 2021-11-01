import numpy

from . import TwoWayDHCross

from pybropt.core.error import check_is_int

class GenerationalTwoWayDHCross(TwoWayDHCross):
    """docstring for GenerationalTwoWayDHCross."""

    def __init__(self, gmult, rng = None, **kwargs):
        """
        Parameters
        ----------
        gmult : int
            Generation multiplier to apply to group labels
        """
        super(GenerationalTwoWayDHCross, self).__init__(
            rng = rng,
            **kwargs
        )
        # check data types
        check_is_int(gmult, "gmult")

        # make assignments
        self.gmult = gmult

    def mate(self, pgmat, sel, ncross, nprogeny, t_cur, s = 0, **kwargs):
        """
        Mate individuals according to a 2-way mate selection scheme.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            A GenotypeMatrix containing candidate breeding individuals.
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
        t_cur : int
            Current time.
        s : int, default = 0
            Number of selfing generations post-cross before double haploids are
            generated.

        Returns
        -------
        progeny : PhasedGenotypeMatrix
            A PhasedGenotypeMatrix of progeny.
        """
        # construct taxa_grp
        taxa_grp = numpy.repeat(
            numpy.repeat(
                numpy.arange(len(sel) // 2, dtype = 'int64'),
                ncross
            ),
            nprogeny
        )
        taxa_grp += t_cur * self.gmult

        # construct progeny
        progeny, misc = super(GenerationalTwoWayDHCross, self).mate(
            pgmat = pgmat,
            sel = sel,
            ncross = ncross,
            nprogeny = nprogeny,
            s = s,
            taxa_grp = taxa_grp,    # add taxa_grp as keyword argument
            **kwargs
        )

        # group progeny taxa
        progeny.group_taxa()

        return progeny, misc
