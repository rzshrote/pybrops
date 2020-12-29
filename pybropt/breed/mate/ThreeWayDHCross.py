import numpy
from .util import mat_mate
from . import MatingOperator
import pybropt.util
import pybropt.popgen.gmat

class ThreeWayDHCross(MatingOperator):
    """docstring for ThreeWayDHCross."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(ThreeWayDHCross, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgvmat, sel, ncross, nprogeny, s = 0):
        """
        Mate individuals according to a 3-way mate selection scheme.

        Parameters
        ----------
        pgvmat : PhasedGenotypeVariantMatrix
            A GenotypeVariantMatrix containing candidate breeding individuals.
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the recurrent parent.
                Second index is the female parent.
                Third index is the male parent.
            Example:
                [1,5,3,8,2,7]
                recurrent = 1,8
                female = 5,2
                male = 3,7
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
        # check data type
        pybropt.popgen.gmat.is_PhasedGenotypeVariantMatrix(pgvmat, "pgvmat")

        # get recurrent, female, and male selections; repeat by ncross
        rsel = numpy.repeat(sel[0::3], ncross)
        fsel = numpy.repeat(sel[1::3], ncross)
        msel = numpy.repeat(sel[2::3], ncross)

        # get pointers to genotypes and crossover probabilities, respectively
        geno = pgvmat.geno
        xoprob = pgvmat.vrnt_xoprob

        # create F1 genotypes
        f1geno = mat_mate(geno, geno, fsel, msel, xoprob)

        # generate selection array for all hybrid lines
        asel = numpy.arange(f1geno.shape[1])

        # generate three way crosses
        hgeno = mat_mate(geno, f1geno, rsel, asel, xoprob)

        # self down hybrids if needed
        for i in range(s):
            # self hybrids
            hgeno = mat_mate(hgeno, hgeno, asel, asel, xoprob)

        # generate doubled haploids
        dhgeno = mat_dh(hgeno, asel, xoprob)

        # create new PhasedGenotypeVariantMatrix
        progeny = PhasedGenotypeVariantMatrix(
            geno = dhgeno,
            vrnt_chrgrp = pgvmat.vrnt_chrgrp,
            vrnt_phypos = pgvmat.vrnt_phypos,
            vrnt_name = pgvmat.vrnt_name,
            vrnt_genpos = pgvmat.vrnt_genpos,
            vrnt_xoprob = pgvmat.vrnt_xoprob,
            vrnt_hapgrp = pgvmat.vrnt_hapgrp,
            vrnt_mask = pgvmat.vrnt_mask
        )

        return progeny

    def mate_cfg(self, mcfg):
        """
        Mate individuals according to a MatingConfiguration.

        Parameters
        ----------
        mcfg : MatingConfiguration
            A MatingConfiguration from which to implement matings.

        Returns
        -------
        progeny : PhasedGenotypeVariantMatrix
            A PhasedGenotypeVariantMatrix of progeny.
        """
        # extract the mating configuration as a dictionary
        d = mcfg.config()

        # unpack and return results
        return self.mate(**d)


################################################################################
################################## Utilities ###################################
################################################################################
def is_ThreeWayDHCross(v):
    return isinstance(v, ThreeWayDHCross)

def check_is_ThreeWayDHCross(v, varname):
    if not isinstance(v, ThreeWayDHCross):
        raise TypeError("'%s' must be a ThreeWayDHCross." % varname)

def cond_check_is_ThreeWayDHCross(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_ThreeWayDHCross(v, varname)
