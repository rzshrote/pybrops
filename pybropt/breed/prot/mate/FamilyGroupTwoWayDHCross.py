import numpy

from . import TwoWayDHCross

from pybropt.core.error import check_is_int

class FamilyGroupTwoWayDHCross(TwoWayDHCross):
    """docstring for FamilyGroupTwoWayDHCross."""

    def __init__(self, counter = 0, rng = None, **kwargs):
        """
        Parameters
        ----------
        counter : int
            Counter used to apply group labels to families of progeny.
        """
        super(FamilyGroupTwoWayDHCross, self).__init__(rng = rng, **kwargs)

        # check data types
        check_is_int(counter, "counter")

        # make assignments
        self.counter = counter

    def mate(self, pgmat, sel, ncross, nprogeny, s = 0, **kwargs):
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
        nfam = len(sel) // 2            # calculate number of families
        stfam = self.counter            # calculate start family numbers
        spfam = self.counter + nfam     # calculate stop family numbers

        # construct taxa_grp
        taxa_grp = numpy.repeat(
            numpy.repeat(
                numpy.arange(stfam, spfam, dtype = 'int64'),
                ncross
            ),
            nprogeny
        )

        # increment counter
        self.counter += nfam

        # construct progeny
        progeny, misc = super(FamilyGroupTwoWayDHCross, self).mate(
            pgmat = pgmat,
            sel = sel,
            ncross = ncross,
            nprogeny = nprogeny,
            s = s,
            taxa_grp = taxa_grp,    # add taxa_grp as keyword argument for object construction
            **kwargs
        )

        # group progeny taxa
        progeny.group_taxa()

        return progeny, misc
