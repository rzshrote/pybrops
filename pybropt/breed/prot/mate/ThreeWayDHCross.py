import numpy
import pybropt.core.random
from pybropt.breed.prot.mate.util import mat_mate
from pybropt.breed.prot.mate.util import mat_dh
from pybropt.breed.prot.mate.MatingProtocol import MatingProtocol
from pybropt.core.error import cond_check_is_Generator
from pybropt.core.error import check_ndarray_len_is_multiple_of_3
from pybropt.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybropt.popgen.gmat.DensePhasedGenotypeMatrix import check_is_DensePhasedGenotypeMatrix

class ThreeWayDHCross(MatingProtocol):
    """docstring for ThreeWayDHCross."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, rng = None, **kwargs):
        """
        Constructor for the concrete class ThreeWayDHCross.

        Parameters
        ----------
        rng : numpy.Generator
            Random number source.
        kwargs : dict
            Additional keyword arguments.
        """
        super(ThreeWayDHCross, self).__init__(**kwargs)

        # check data types
        cond_check_is_Generator(rng, "rng")

        # make assignments
        self.rng = pybropt.core.random if rng is None else rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgmat, sel, ncross, nprogeny, miscout = None, s = 0, **kwargs):
        """
        Mate individuals according to a 3-way mate selection scheme.

        Parameters
        ----------
        pgmat : DensePhasedGenotypeMatrix
            A GenotypeMatrix containing candidate breeding individuals.
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
        check_ndarray_len_is_multiple_of_3(sel, "sel")

        # get recurrent, female, and male selections; repeat by ncross
        rsel = numpy.repeat(sel[0::3], ncross)  # recurrent parent
        fsel = numpy.repeat(sel[1::3], ncross)  # female parent
        msel = numpy.repeat(sel[2::3], ncross)  # male parent

        # get pointers to genotypes and crossover probabilities, respectively
        geno = pgmat.mat
        xoprob = pgmat.vrnt_xoprob

        # create F1 genotypes
        f1geno = mat_mate(geno, geno, fsel, msel, xoprob, self.rng)

        # generate selection array for all hybrid lines
        hsel = numpy.arange(f1geno.shape[1])

        # create backcross genotypes
        bcgeno = mat_mate(geno, f1geno, rsel, hsel, xoprob, self.rng)

        # generate selection array for all backcross lines
        bcsel = numpy.arange(bcgeno.shape[1])

        # self down hybrids if needed
        for i in range(s):
            # self backcross lines
            bcgeno = mat_mate(bcgeno, bcgeno, bcsel, bcsel, xoprob, self.rng)

        # generate selection array for progeny
        psel = numpy.repeat(bcsel, nprogeny)

        # generate doubled haploids
        dhgeno = mat_dh(bcgeno, psel, xoprob, self.rng)

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
            kwargs : dict
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
def is_ThreeWayDHCross(v):
    """
    Determine whether an object is a ThreeWayDHCross.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a ThreeWayDHCross object instance.
    """
    return isinstance(v, ThreeWayDHCross)

def check_is_ThreeWayDHCross(v, varname):
    """
    Check if object is of type ThreeWayDHCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, ThreeWayDHCross):
        raise TypeError("'%s' must be a ThreeWayDHCross." % varname)

def cond_check_is_ThreeWayDHCross(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type ThreeWayDHCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a ThreeWayDHCross.
    """
    if cond(v):
        check_is_ThreeWayDHCross(v, varname)
