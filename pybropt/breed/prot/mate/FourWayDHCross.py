import numpy
import pybropt.core.random
from pybropt.breed.prot.mate.util import mat_mate
from pybropt.breed.prot.mate.util import mat_dh
from pybropt.breed.prot.mate.MatingProtocol import MatingProtocol
from pybropt.core.error import cond_check_is_Generator
from pybropt.core.error import check_ndarray_len_is_multiple_of_4
from pybropt.popgen.gmat import check_is_DensePhasedGenotypeMatrix

class FourWayDHCross(MatingProtocol):
    """docstring for FourWayDHCross."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, rng = None, **kwargs):
        """
        Constructor for the concrete class FourWayDHCross.

        Parameters
        ----------
        rng : numpy.Generator
            Random number source.
        kwargs : dict
            Additional keyword arguments.
        """
        super(FourWayDHCross, self).__init__(**kwargs)

        # check data types
        cond_check_is_Generator(rng, "rng")

        # make assignments
        self.rng = pybropt.core.random if rng is None else rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgmat, sel, ncross, nprogeny, miscout = None, s = 0, **kwargs):
        """
        Mate individuals according to a 4-way mate selection scheme.

        Parameters
        ----------
        pgmat : DensePhasedGenotypeMatrix
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
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        s : int, default = 0
            Number of selfing generations post-cross before double haploids are
            generated.
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
        check_ndarray_len_is_multiple_of_4(sel, "sel")

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

        # generate selection array for hybrid lines
        absel = numpy.arange(abgeno.shape[1])
        cdsel = numpy.arange(cdgeno.shape[1])

        # generate dihybrid cross
        dihgeno = mat_mate(abgeno, cdgeno, absel, cdsel, xoprob, self.rng)

        # generate selection array for dihybrid lines
        dihsel = numpy.arange(dihgeno.shape[1])

        # self down hybrids if needed
        for i in range(s):
            # self hybrids
            dihgeno = mat_mate(dihgeno, dihgeno, dihsel, dihsel, xoprob, self.rng)

        # generate selection array for progeny
        psel = numpy.repeat(dihsel, nprogeny)

        # generate doubled haploids
        dhgeno = mat_dh(dihgeno, psel, xoprob, self.rng)

        # create new DensePhasedGenotypeMatrix
        progeny = pgmat.__class__(
            mat = dhgeno,
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
def is_FourWayDHCross(v):
    """
    Determine whether an object is a FourWayDHCross.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a FourWayDHCross object instance.
    """
    return isinstance(v, FourWayDHCross)

def check_is_FourWayDHCross(v, varname):
    """
    Check if object is of type FourWayDHCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, FourWayDHCross):
        raise TypeError("'%s' must be a FourWayDHCross." % varname)

def cond_check_is_FourWayDHCross(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type FourWayDHCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a FourWayDHCross.
    """
    if cond(v):
        check_is_FourWayDHCross(v, varname)
