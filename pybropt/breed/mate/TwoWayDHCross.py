import numpy
from . import mat_mate
from . import mat_dh
from . import MatingOperator
from pybropt.popgen.gmat import DensePhasedGenotypeVariantMatrix
from pybropt.popgen.gmat import check_is_DensePhasedGenotypeVariantMatrix

class TwoWayDHCross(MatingOperator):
    """docstring for TwoWayDHCross."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, rng, **kwargs):
        super(TwoWayDHCross, self).__init__(**kwargs)
        self.rng = rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, t_cur, t_max, pgvmat, sel, ncross, nprogeny, s = 0, **kwargs):
        """
        Mate individuals according to a 2-way mate selection scheme.

        Parameters
        ----------
        pgvmat : DensePhasedGenotypeVariantMatrix
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
        progeny : DensePhasedGenotypeVariantMatrix
            A DensePhasedGenotypeVariantMatrix of progeny.
        """
        # print(id(pgvmat))
        # check data type
        check_is_DensePhasedGenotypeVariantMatrix(pgvmat, "pgvmat")

        # get female and male selections; repeat by ncross
        fsel = numpy.repeat(sel[0::2], ncross)
        msel = numpy.repeat(sel[1::2], ncross)

        # get pointers to genotypes and crossover probabilities, respectively
        geno = pgvmat.mat
        xoprob = pgvmat.vrnt_xoprob

        # create hybrid genotypes
        hgeno = mat_mate(geno, geno, fsel, msel, xoprob, self.rng)

        # generate selection array for all hybrid lines
        asel = numpy.repeat(numpy.arange(hgeno.shape[1]), nprogeny)

        # self down hybrids if needed
        for i in range(s):
            # self hybrids
            hgeno = mat_mate(hgeno, hgeno, asel, asel, xoprob, self.rng)

        # generate doubled haploids
        dhgeno = mat_dh(hgeno, asel, xoprob, self.rng)

        # create new DensePhasedGenotypeVariantMatrix
        progeny = DensePhasedGenotypeVariantMatrix(
            mat = dhgeno,
            vrnt_chrgrp = pgvmat.vrnt_chrgrp,
            vrnt_phypos = pgvmat.vrnt_phypos,
            vrnt_name = pgvmat.vrnt_name,
            vrnt_genpos = pgvmat.vrnt_genpos,
            vrnt_xoprob = pgvmat.vrnt_xoprob,
            vrnt_hapgrp = pgvmat.vrnt_hapgrp,
            vrnt_mask = pgvmat.vrnt_mask
        )

        progeny.vrnt_chrgrp_name = pgvmat.vrnt_chrgrp_name
        progeny.vrnt_chrgrp_stix = pgvmat.vrnt_chrgrp_stix
        progeny.vrnt_chrgrp_spix = pgvmat.vrnt_chrgrp_spix
        progeny.vrnt_chrgrp_len = pgvmat.vrnt_chrgrp_len

        misc = {}

        return progeny, misc



################################################################################
################################## Utilities ###################################
################################################################################
def is_TwoWayDHCross(v):
    return isinstance(v, TwoWayDHCross)

def check_is_TwoWayDHCross(v, varname):
    if not isinstance(v, TwoWayDHCross):
        raise TypeError("'%s' must be a TwoWayDHCross." % varname)

def cond_check_is_TwoWayDHCross(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_TwoWayDHCross(v, varname)
