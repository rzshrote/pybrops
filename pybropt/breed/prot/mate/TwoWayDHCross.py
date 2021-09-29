import numpy
import pybropt.core.random
from . import mat_mate
from . import mat_dh
from . import MatingProtocol
from pybropt.core.error import cond_check_is_Generator
from pybropt.core.error import check_ndarray_len_is_multiple_of_2
from pybropt.popgen.gmat import DensePhasedGenotypeMatrix
from pybropt.popgen.gmat import check_is_DensePhasedGenotypeMatrix

class TwoWayDHCross(MatingProtocol):
    """docstring for TwoWayDHCross."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, rng = None, **kwargs):
        """
        Constructor for the concrete class TwoWayDHCross.

        Parameters
        ----------
        rng : numpy.Generator
            Random number source.
        **kwargs : dict
            Additional keyword arguments.
        """
        super(TwoWayDHCross, self).__init__(**kwargs)

        # check data types
        cond_check_is_Generator(rng, "rng")

        # make assignments
        self.rng = pybropt.core.random if rng is None else rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgmat, sel, ncross, nprogeny, s = 0, **kwargs):
        """
        Mate individuals according to a 2-way mate selection scheme, then create
        doubled haploid progenies.

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
        check_ndarray_len_is_multiple_of_2(sel, "sel")

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

        # generate doubled haploids
        dhgeno = mat_dh(hgeno, asel, xoprob, self.rng)

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

        # empty additional output
        misc = {}

        return progeny, misc



################################################################################
################################## Utilities ###################################
################################################################################
def is_TwoWayDHCross(v):
    """
    Determine whether an object is a TwoWayDHCross.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a TwoWayDHCross object instance.
    """
    return isinstance(v, TwoWayDHCross)

def check_is_TwoWayDHCross(v, varname):
    """
    Check if object is of type TwoWayDHCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, TwoWayDHCross):
        raise TypeError("'%s' must be a TwoWayDHCross." % varname)

def cond_check_is_TwoWayDHCross(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type TwoWayDHCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a TwoWayDHCross.
    """
    if cond(v):
        check_is_TwoWayDHCross(v, varname)
