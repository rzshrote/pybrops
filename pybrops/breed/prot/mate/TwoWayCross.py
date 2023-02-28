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
from pybrops.core.error.error_type_numpy import check_is_Integral_or_ndarray, check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral, check_is_dict, check_is_int
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
            progeny_counter: Integral = 0, 
            family_counter: Integral = 0, 
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
    def progeny_counter(self) -> Integral:
        """Description for property progeny_counter."""
        return self._progeny_counter
    @progeny_counter.setter
    def progeny_counter(self, value: Integral) -> None:
        """Set data for property progeny_counter."""
        check_is_Integral(value, "progeny_counter")
        self._progeny_counter = value
    @progeny_counter.deleter
    def progeny_counter(self) -> None:
        """Delete data for property progeny_counter."""
        del self._progeny_counter

    @property
    def family_counter(self) -> Integral:
        """Description for property family_counter."""
        return self._family_counter
    @family_counter.setter
    def family_counter(self, value: Integral) -> None:
        """Set data for property family_counter."""
        check_is_Integral(value, "family_counter")
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
            ncross: Union[Integral,numpy.ndarray], 
            nprogeny: Union[Integral,numpy.ndarray], 
            miscout: Optional[dict] = None, 
            nself: Integral = 0, 
            **kwargs: dict
        ) -> PhasedGenotypeMatrix:
        """
        Mate individuals according to a 2-way mate selection scheme.

        Example crossing diagram::

                                pgmat
                                │                       sel = [A,B,...]
                               AxB
                    ┌───────────┴───────────┐           ncross * nprogeny = 2
                   AxB                     AxB          duplicate cross (ncross * nprogeny) times
                    │                       │           nself = 2
                S0(AxB)                 S0(AxB)         first self
                    │                       │
                S1(AxB)                 S1(AxB)         second self
                    |                       |
                S1(AxB)                 S1(AxB)         final result

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
            Number of progeny to generate per cross.
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
        # check data types
        check_is_DensePhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_ndarray(sel, "sel")
        check_ndarray_len_is_multiple_of_2(sel, "sel")
        check_is_Integral_or_ndarray(ncross, "ncross")
        check_is_Integral_or_ndarray(nprogeny, "nprogeny")
        if miscout is not None:
            check_is_dict(miscout, "miscout")
        check_is_Integral(nself, "nself")

        # get female and male selections; repeat by ncross
        fsel = numpy.repeat(sel[0::2], ncross * nprogeny)
        msel = numpy.repeat(sel[1::2], ncross * nprogeny)

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

        ########################################################################
        ######################### Metadata generation ##########################
        # generate line names
        progcnt = hgeno.shape[1]                # get number of hybrid progeny generated
        riter = range(                          # range iterator for line names
            self.progeny_counter,               # start progeny number (inclusive)
            self.progeny_counter + progcnt      # stop progeny number (exclusive)
        )
        # create taxa names
        taxa = numpy.array(["2w"+str(i).zfill(7) for i in riter], dtype = "object")
        self.progeny_counter += progcnt         # increment counter

        # calculate taxa family groupings
        nfam = len(sel) // 2                    # calculate number of families
        taxa_grp = numpy.repeat(                # construct taxa_grp vector
            numpy.arange(                       # repeat for crosses * progeny
                self.family_counter,            # start family number (inclusive)
                self.family_counter + nfam,     # stop family number (exclusive)
                dtype = 'int64'
            ),
            ncross * nprogeny
        )
        self.family_counter += nfam             # increment counter
        
        ########################################################################
        ########################## Output generation ###########################
        # create new DensePhasedGenotypeMatrix
        progeny = pgmat.__class__(
            mat = hgeno,
            taxa = taxa,
            taxa_grp = taxa_grp,
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

        # group progeny taxa
        progeny.group_taxa()

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
