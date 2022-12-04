"""
Module implementing mating protocols for four-way crosses.
"""

import numpy

from pybrops.breed.prot.mate.util import mat_dh
from pybrops.breed.prot.mate.util import mat_mate
from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import check_is_DensePhasedGenotypeMatrix
from pybrops.core.random.prng import global_prng

class FourWayCross(MatingProtocol):
    """
    Class implementing mating protocols for four-way crosses.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, rng = None, **kwargs):
        super(FourWayCross, self).__init__(**kwargs)

        # make assignments
        self.rng = rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def rng():
        doc = "The rng property."
        def fget(self):
            """Get value for rng."""
            return self._rng
        def fset(self, value):
            """Set value for rng."""
            if value is None:
                check_is_Generator_or_RandomState(value, "rng")
            else:
                value = global_prng
            self._rng = value
        def fdel(self):
            """Delete value for rng."""
            del self._rng
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    rng = property(**rng())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgmat: DensePhasedGenotypeMatrix, sel: numpy.ndarray, ncross: int, nprogeny: int, miscout = None, s = 0, **kwargs):
        """
        Mate individuals according to a 4-way mate selection scheme.

        Parameters
        ----------
        pgmat : DensePhasedGenotypeMatrix
            A GenotypeMatrix containing candidate breeding individuals.
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape ``(k,)``.

            Where:

            - ``k`` is the number of selected individuals.

            Indices are paired as follows:

            - First index is the female parent 2.
            - Second index is the male parent 2.
            - Third index is the female parent 1.
            - Fourth index is the male parent 1.

            Example::

                sel = [1,5,3,8]
                female2 = 1
                male2 = 5
                female1 = 3
                male1 = 8
        ncross : numpy.ndarray
            Number of cross patterns to perform.
        nprogeny : numpy.ndarray
            Number of doubled haploid progeny to generate per cross.
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        s : int, default = 0
            Number of selfing generations post-cross.
        kwargs : dict
            Additional keyword arguments to be passed to constructor for the
            output DensePhasedGenotypeMatrix.

        Returns
        -------
        out : PhasedGenotypeMatrix
            A PhasedGenotypeMatrix of progeny.
        """
        # check data type
        check_is_DensePhasedGenotypeMatrix(pgmat, "pgmat")

        # get female2, male2, female1, and male1 selections; repeat by ncross
        f2sel = numpy.repeat(sel[0::4], ncross)
        m2sel = numpy.repeat(sel[1::4], ncross)
        f1sel = numpy.repeat(sel[2::4], ncross)
        m1sel = numpy.repeat(sel[3::4], ncross)

        # get pointers to genotypes and crossover probabilities, respectively
        geno = pgmat.mat
        xoprob = pgmat.vrnt_xoprob

        # create F1 genotypes
        abgeno = mat_mate(geno, geno, f1sel, m1sel, xoprob, self.rng)
        cdgeno = mat_mate(geno, geno, f2sel, m2sel, xoprob, self.rng)

        # generate selection array for all hybrid lines
        asel = numpy.repeat(numpy.arange(abgeno.shape[1]), nprogeny)

        # generate dihybrid cross
        hgeno = mat_mate(abgeno, cdgeno, asel, asel, xoprob, self.rng)

        # self down hybrids if needed
        for i in range(s):
            # self hybrids
            hgeno = mat_mate(hgeno, hgeno, asel, asel, xoprob, self.rng)

        # create new DensePhasedGenotypeMatrix
        progeny = pgmat.__class__(
            mat = hgeno,
            vrnt_chrgrp = pgmat.vrnt_chrgrp,
            vrnt_phypos = pgmat.vrnt_phypos,
            vrnt_name = pgmat.vrnt_name,
            vrnt_genpos = pgmat.vrnt_genpos,
            vrnt_xoprob = pgmat.vrnt_xoprob,
            vrnt_hapgrp = pgmat.vrnt_hapgrp,
            vrnt_mask = pgmat.vrnt_mask,
            **kwargs
        )

        # copy metadata
        progeny.vrnt_chrgrp_name = pgmat.vrnt_chrgrp_name
        progeny.vrnt_chrgrp_stix = pgmat.vrnt_chrgrp_stix
        progeny.vrnt_chrgrp_spix = pgmat.vrnt_chrgrp_spix
        progeny.vrnt_chrgrp_len = pgmat.vrnt_chrgrp_len

        return progeny



################################################################################
################################## Utilities ###################################
################################################################################
def is_FourWayCross(v):
    return isinstance(v, FourWayCross)

def check_is_FourWayCross(v, varname):
    if not isinstance(v, FourWayCross):
        raise TypeError("'%s' must be a FourWayCross." % varname)

def cond_check_is_FourWayCross(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_FourWayCross(v, varname)
