import numpy
import pybrops.core.random
from pybrops.breed.prot.mate.util import mat_mate
from pybrops.breed.prot.mate.util import mat_dh
from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol
from pybrops.core.error import cond_check_is_Generator
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import check_is_DensePhasedGenotypeMatrix

class SelfCross(MatingProtocol):
    """docstring for SelfCross."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, rng = None, **kwargs):
        """
        Constructor for the concrete class SelfCross.

        Parameters
        ----------
        rng : numpy.Generator
            Random number source.
        kwargs : dict
            Additional keyword arguments.
        """
        super(SelfCross, self).__init__(**kwargs)

        # check data types
        cond_check_is_Generator(rng, "rng")

        # make assignments
        self.rng = pybrops.core.random if rng is None else rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgmat, sel, ncross, nprogeny, miscout = None, s = 0, **kwargs):
        """
        Self-fertilize individuals.

        Example crossing diagram::

                     pgmat
                       │                sel = [A,...]
                       A
                 ┌─────┴─────┐          ncross = 2
                 A           A
              ┌──┴──┐     ┌──┴──┐       initial self, nprogeny = 2
            S0(A) S0(A) S0(A) S0(A)     cross pattern finished.
              │     │     │     │       s = 2
            S1(A) S1(A) S1(A) S1(A)     first self
              │     │     │     │
            S2(A) S2(A) S2(A) S2(A)     second self, final result

        Parameters
        ----------
        pgmat : DensePhasedGenotypeMatrix
            A DensePhasedGenotypeMatrix containing candidate breeding
            individuals.
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape ``(k,)``.

            Where:

            - ``k`` is the number of selected individuals.

            Indices are as follows:

            - All indices are female.

            Example::

                [1,5,3,8,2,7]
                female = 1,5,3,8,2,7
        ncross : numpy.ndarray
            Number of cross patterns to perform.
        nprogeny : numpy.ndarray
            Number of progeny to generate per cross.
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        s : int, default = 0
            Number of generations of single seed descent post-cross pattern.
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

        # get pointers to genotypes and crossover probabilities, respectively
        geno = pgmat.mat
        xoprob = pgmat.vrnt_xoprob

        # get female selections; repeat by ncross
        fsel = numpy.repeat(sel, ncross*nprogeny)

        # self genotypes
        sgeno = mat_mate(geno, geno, fsel, fsel, xoprob, self.rng)

        # perform single seed descent
        for i in range(s):
            # generate selection array for all selfed lines
            ssel = numpy.arange(sgeno.shape[1])
            # self lines
            sgeno = mat_mate(sgeno, sgeno, ssel, ssel, xoprob, self.rng)

        # create new DensePhasedGenotypeMatrix
        progeny = pgmat.__class__(
            mat = sgeno,
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
def is_SelfCross(v):
    """
    Determine whether an object is a SelfCross.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a SelfCross object instance.
    """
    return isinstance(v, SelfCross)

def check_is_SelfCross(v, varname):
    """
    Check if object is of type SelfCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SelfCross):
        raise TypeError("'%s' must be a SelfCross." % varname)

def cond_check_is_SelfCross(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type SelfCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a SelfCross.
    """
    if cond(v):
        check_is_SelfCross(v, varname)
