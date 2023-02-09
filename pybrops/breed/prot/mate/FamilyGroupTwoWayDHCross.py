"""
Module implementing mating protocols for two-way DH crosses and storing family
information.
"""

import numpy

from pybrops.breed.prot.mate.TwoWayDHCross import TwoWayDHCross
from pybrops.core.error import check_is_int

class FamilyGroupTwoWayDHCross(TwoWayDHCross):
    """
    Class implementing mating protocols for two-way DH crosses while storing
    family information in the resulting genotype matrix.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, counter = 0, rng = None, **kwargs: dict):
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

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgmat, sel, ncross, nprogeny, miscout = None, s = 0, **kwargs: dict):
        """
        Mate individuals according to a 2-way mate selection scheme.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            A GenotypeMatrix containing candidate breeding individuals.
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape ``(k,)``.

            Where:

            - ``k`` is the number of selected individuals.

            Indices are paired as follows:

            - Even indices are female.
            - Odd indices are male.

            Example::

                sel = [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7
        ncross : numpy.ndarray
            Number of cross patterns to perform.
        nprogeny : numpy.ndarray
            Number of doubled haploid progeny to generate per cross.
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        s : int, default = 0
            Number of selfing generations post-cross before double haploids are
            generated.

        Returns
        -------
        out : PhasedGenotypeMatrix
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
        progeny = super(FamilyGroupTwoWayDHCross, self).mate(
            pgmat = pgmat,
            sel = sel,
            ncross = ncross,
            nprogeny = nprogeny,
            miscout = miscout,
            s = s,
            taxa_grp = taxa_grp,    # add taxa_grp as keyword argument for object construction
            **kwargs
        )

        # group progeny taxa
        progeny.group_taxa()

        return progeny
