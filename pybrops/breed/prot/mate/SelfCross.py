"""
Module implementing mating protocols for self-fertilization.
"""

from typing import Any
import numpy

from pybrops.breed.prot.mate.util import mat_mate
from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import check_is_DensePhasedGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.core.random.prng import global_prng

class SelfCross(MatingProtocol):
    """
    Class implementing mating protocols for self-fertilization.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, rng = None, **kwargs):
        """
        Constructor for the concrete class SelfCross.

        Parameters
        ----------
        rng : numpy.random.Generator, numpy.random.RandomState
            Random number source.
        kwargs : dict
            Additional keyword arguments.
        """
        super(SelfCross, self).__init__(**kwargs)

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
                value = global_prng
            check_is_Generator_or_RandomState(value, "rng")
            self._rng = value
        def fdel(self):
            """Delete value for rng."""
            del self._rng
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    rng = property(**rng())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgmat: DensePhasedGenotypeMatrix, sel: numpy.ndarray, ncross: int, nprogeny: int, miscout: dict = None, s: int = 0, **kwargs):
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
def is_SelfCross(v: Any) -> bool:
    """
    Determine whether an object is a SelfCross.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a SelfCross object instance.
    """
    return isinstance(v, SelfCross)

def check_is_SelfCross(v: Any, varname: str) -> None:
    """
    Check if object is of type SelfCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SelfCross):
        raise TypeError("'%s' must be a SelfCross." % varname)
