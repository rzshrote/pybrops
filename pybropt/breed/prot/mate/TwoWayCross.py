import numpy
import pybropt.core.random
from . import mat_mate
from . import mat_dh
from . import MatingProtocol
from pybropt.core.error import cond_check_is_Generator
from pybropt.popgen.gmat import DensePhasedGenotypeMatrix
from pybropt.popgen.gmat import check_is_DensePhasedGenotypeMatrix

class TwoWayCross(MatingProtocol):
    """docstring for TwoWayCross."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, rng = None, **kwargs):
        super(TwoWayCross, self).__init__(**kwargs)

        # check data types
        cond_check_is_Generator(rng, "rng")

        # make assignments
        self.rng = pybropt.core.random if rng is None else rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgmat, sel, ncross, nprogeny, s = 0, **kwargs):
        """
        Mate individuals according to a 2-way mate selection scheme.

        Parameters
        ----------
        pgmat : DensePhasedGenotypeMatrix
            A DensePhasedGenotypeMatrix containing candidate breeding
            individuals.
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
        **kwargs : dict
            Additional keyword arguments to be passed to constructor for the
            output DensePhasedGenotypeMatrix.

        Returns
        -------
        progeny : DensePhasedGenotypeMatrix
            A DensePhasedGenotypeMatrix of progeny. Output is of the same class
            as 'pgmat', but are guaranteed to be a DensePhasedGenotypeMatrix.
        """
        # check data type
        check_is_DensePhasedGenotypeMatrix(pgmat, "pgmat")

        # get female and male selections; repeat by ncross
        fsel = numpy.repeat(sel[0::2], ncross)
        msel = numpy.repeat(sel[1::2], ncross)

        # get pointers to genotypes and crossover probabilities, respectively
        geno = pgmat.mat
        xoprob = pgmat.vrnt_xoprob

        # create hybrid genotypes
        hgeno = mat_mate(geno, geno, fsel, msel, xoprob, self.rng)

        # generate selection array for all hybrid lines
        asel = numpy.repeat(numpy.arange(hgeno.shape[1]), nprogeny)

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

        # empty additional output
        misc = {}

        return progeny, misc



################################################################################
################################## Utilities ###################################
################################################################################
def is_TwoWayCross(v):
    return isinstance(v, TwoWayCross)

def check_is_TwoWayCross(v, varname):
    if not isinstance(v, TwoWayCross):
        raise TypeError("'%s' must be a TwoWayCross." % varname)

def cond_check_is_TwoWayCross(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_TwoWayCross(v, varname)