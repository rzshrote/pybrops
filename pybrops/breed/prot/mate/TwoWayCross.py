"""
Module implementing mating protocols for two-way crosses.
"""

from numbers import Integral
from typing import Any, Optional, Union
import numpy

from pybrops.breed.prot.mate.util import mat_mate
from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol
from pybrops.core.error import check_ndarray_len_is_multiple_of_2
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_type_python import check_is_int
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import check_is_DensePhasedGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.core.random.prng import global_prng
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class TwoWayCross(MatingProtocol):
    """
    Class implementing mating protocols for two-way crosses.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            progeny_counter: int = 0, 
            family_counter: int = 0, 
            rng: Union[numpy.random.Generator,numpy.random.RandomState,None] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class TwoWayCross.

        Parameters
        ----------
        rng : numpy.random.Generator, numpy.random.RandomState
            Random number source.
        kwargs : dict
            Additional keyword arguments.
        """
        super(TwoWayCross, self).__init__(**kwargs)

        # make assignments
        self.progeny_counter = progeny_counter
        self.family_counter = family_counter
        self.rng = rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def nparent(self) -> Integral:
        """Number of parents the mating protocol requires."""
        return 2
    @nparent.setter
    def nparent(self, value: Integral) -> None:
        """Set number of parents the mating protocol requires."""
        error_readonly("nparent")
    @nparent.deleter
    def nparent(self) -> None:
        """Delete number of parents the mating protocol requires."""
        error_readonly("nparent")

    @property
    def progeny_counter(self) -> int:
        """Description for property progeny_counter."""
        return self._progeny_counter
    @progeny_counter.setter
    def progeny_counter(self, value: int) -> None:
        """Set data for property progeny_counter."""
        check_is_int(value, "progeny_counter")
        self._progeny_counter = value
    @progeny_counter.deleter
    def progeny_counter(self) -> None:
        """Delete data for property progeny_counter."""
        del self._progeny_counter

    @property
    def family_counter(self) -> int:
        """Description for property family_counter."""
        return self._family_counter
    @family_counter.setter
    def family_counter(self, value: int) -> None:
        """Set data for property family_counter."""
        check_is_int(value, "family_counter")
        self._family_counter = value
    @family_counter.deleter
    def family_counter(self) -> None:
        """Delete data for property family_counter."""
        del self._family_counter

    @property
    def rng(self) -> Union[numpy.random.Generator,numpy.random.RandomState]:
        """Random number generator."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[numpy.random.Generator,numpy.random.RandomState,None]) -> None:
        """Set random number generator."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng")
        self._rng = value
    @rng.deleter
    def rng(self) -> None:
        """Delete random number generator."""
        del self._rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            sel: numpy.ndarray, 
            ncross: Union[int,numpy.ndarray], 
            nprogeny: Union[int,numpy.ndarray], 
            miscout: Optional[dict] = None, 
            nself: int = 0, 
            **kwargs: dict
        ) -> PhasedGenotypeMatrix:
        """
        Mate individuals according to a 2-way mate selection scheme.

        Parameters
        ----------
        pgmat : DensePhasedGenotypeMatrix
            A DensePhasedGenotypeMatrix containing candidate breeding
            individuals.
        sel : numpy.ndarray
            A 1D array of indices of selected individuals of shape ``(k,)``.

            Where:

            - ``k`` is the number of selected individuals.

            Indices are paired as follows:

            - Even indices are female.
            - Odd indices are male.

            Example::

                sel = [1,5,3,8,2,7]
                female = 1,3,2
                male = 5,8,7
        ncross : numpy.ndarray
            Number of cross patterns to perform.
        nprogeny : numpy.ndarray
            Number of doubled haploid progeny to generate per cross.
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
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
        for i in range(nself):
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
def is_TwoWayCross(v: Any) -> bool:
    """
    Determine whether an object is a TwoWayCross.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a TwoWayCross object instance.
    """
    return isinstance(v, TwoWayCross)

def check_is_TwoWayCross(v: Any, varname: str) -> None:
    """
    Check if object is of type TwoWayCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, TwoWayCross):
        raise TypeError("'%s' must be a TwoWayCross." % varname)
