import numpy
from . import MatingOperator
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import check_is_DensePhasedGenotypeMatrix

class FourWayDHCross(MatingOperator):
    """docstring for FourWayDHCross."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict):
        super(FourWayDHCross, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgvmat, sel, ncross, nprogeny, s = 0):
        """
        Mate individuals according to a 4-way mate selection scheme.

        Parameters
        ----------
        pgvmat : DensePhasedGenotypeMatrix
            A GenotypeMatrix containing candidate breeding individuals.
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape (k,).
            Where:
                'k' is the number of selected individuals.
            Indices are paired as follows:
                First index is the female parent 2.
                Second index is the male parent 2.
                Third index is the female parent 1.
                Fourth index is the male parent 1.
            Example:
                [1,5,3,8]
                female2 = 1
                male2 = 5
                female1 = 3
                male1 = 8
        ncross : numpy.ndarray
            Number of cross patterns to perform.
        nprogeny : numpy.ndarray
            Number of doubled haploid progeny to generate per cross.
        s : int, default = 0
            Number of selfing generations post-cross before double haploids are
            generated.

        Returns
        -------
        progeny : DensePhasedGenotypeMatrix
            A DensePhasedGenotypeMatrix of progeny.
        """
        # check data type
        check_is_DensePhasedGenotypeMatrix(pgvmat, "pgvmat")

        # get female2, male2, female1, and male1 selections; repeat by ncross
        f2sel = numpy.repeat(sel[0::4], ncross)
        m2sel = numpy.repeat(sel[1::4], ncross)
        f1sel = numpy.repeat(sel[2::4], ncross)
        m1sel = numpy.repeat(sel[3::4], ncross)

        # get pointers to genotypes and crossover probabilities, respectively
        geno = pgvmat.geno
        xoprob = pgvmat.vrnt_xoprob

        # create F1 genotypes
        abgeno = mat_mate(geno, geno, f1sel, m1sel, xoprob)
        cdgeno = mat_mate(geno, geno, f2sel, m2sel, xoprob)

        # generate selection array for all hybrid lines
        asel = numpy.repeat(numpy.arange(abgeno.shape[1]), nprogeny)

        # generate dihybrid cross
        hgeno = mat_mate(abgeno, cdgeno, asel, asel, xoprob)

        # self down hybrids if needed
        for i in range(s):
            # self hybrids
            hgeno = mat_mate(hgeno, hgeno, asel, asel, xoprob)

        # generate doubled haploids
        dhgeno = mat_dh(hgeno, asel, xoprob)

        # create new DensePhasedGenotypeMatrix
        progeny = DensePhasedGenotypeMatrix(
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
        progeny : DensePhasedGenotypeMatrix
            A DensePhasedGenotypeMatrix of progeny.
        """
        # extract the mating configuration as a dictionary
        d = mcfg.config()

        # unpack and return results
        return self.mate(**d)



################################################################################
################################## Utilities ###################################
################################################################################
def is_FourWayDHCross(v):
    return isinstance(v, FourWayDHCross)

def check_is_FourWayDHCross(v, varname):
    if not isinstance(v, FourWayDHCross):
        raise TypeError("'%s' must be a FourWayDHCross." % varname)

def cond_check_is_FourWayDHCross(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_FourWayDHCross(v, varname)
