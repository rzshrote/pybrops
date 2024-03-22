"""
Module implementing mating protocols for two-way crosses.
"""

from numbers import Integral
from typing import Optional
from typing import Union
import numpy
from numpy.random import Generator
from numpy.random import RandomState
from pybrops.breed.prot.mate.util import mat_mate
from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_numpy import check_ndarray_shape_eq
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import check_DensePhasedGenotypeMatrix_has_vrnt_xoprob
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import check_is_DensePhasedGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.core.random.prng import global_prng

class TwoWayCross(
        MatingProtocol,
    ):
    """
    Class implementing mating protocols for two-way crosses.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            progeny_counter: Integral = 0, 
            family_counter: Integral = 0, 
            rng: Optional[Union[Generator,RandomState]] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class TwoWayCross.

        Parameters
        ----------
        progeny_counter : Integral
            Progeny counter. This helps create progeny names.
        family_counter : Integral
            Family counter. This helps label groups of progenies as originating 
            from the same family.
        rng : numpy.random.Generator, numpy.random.RandomState
            Random number source.
        kwargs : dict
            Additional keyword arguments.
        """
        # make assignments
        self.progeny_counter = progeny_counter
        self.family_counter = family_counter
        self.rng = rng

    ############################ Object Properties #############################
    @property
    def nparent(self) -> Integral:
        """Number of parents the mating protocol requires."""
        return 2
    @nparent.setter
    def nparent(self, value: Integral) -> None:
        """Set number of parents the mating protocol requires."""
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

    @property
    def family_counter(self) -> Integral:
        """Description for property family_counter."""
        return self._family_counter
    @family_counter.setter
    def family_counter(self, value: Integral) -> None:
        """Set data for property family_counter."""
        check_is_Integral(value, "family_counter")
        self._family_counter = value

    @property
    def rng(self) -> Union[Generator,RandomState]:
        """Random number generator."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[Generator,RandomState,None]) -> None:
        """Set random number generator."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng")
        self._rng = value

    ############################## Object Methods ##############################
    def mate(
            self, 
            pgmat: DensePhasedGenotypeMatrix, 
            xconfig: numpy.ndarray, 
            nmating: Union[Integral,numpy.ndarray], 
            nprogeny: Union[Integral,numpy.ndarray], 
            miscout: Optional[dict] = None, 
            nself: Integral = 0, 
            **kwargs: dict
        ) -> DensePhasedGenotypeMatrix:
        """
        Mate individuals according to a 2-way mate selection scheme.

        Example crossing diagram::

                                pgmat
                                │                       sel = [A,B,...]
                               AxB
                    ┌───────────┴───────────┐           nmating * nprogeny = 2
                   AxB                     AxB          duplicate cross (nmating * nprogeny) times
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
        xconfig : numpy.ndarray
            Array of shape ``(ncross,nparent)`` containing indices specifying a cross
            configuration. Each index corresponds to an individual in ``pgmat``.

            Where:

            - ``ncross`` is the number of crosses to perform.
            - ``nparent`` is the number of parents required for a cross.

            Indices are paired as follows:

            - Even indices are female.
            - Odd indices are male.

            Example::

                xconfig = [[ 1, 5 ],
                           [ 3, 8 ],
                           [ 2, 7 ],
                           ...,
                           [ F, M ]]
                female = [1, 3, 2, ..., F]
                male = [5, 8, 7, ..., M]
        ncross : numpy.ndarray
            Number of cross patterns to perform.
        nprogeny : numpy.ndarray
            Number of progeny to generate per cross.
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        nself : int, default = 0
            Number of selfing generations post-cross before progeny lines are
            generated.
        kwargs : dict
            Additional keyword arguments to be passed to constructor for the
            output DensePhasedGenotypeMatrix.

        Returns
        -------
        out : DensePhasedGenotypeMatrix
            A DensePhasedGenotypeMatrix of progeny.
        """
        # check pgmat
        check_is_DensePhasedGenotypeMatrix(pgmat, "pgmat")
        check_DensePhasedGenotypeMatrix_has_vrnt_xoprob(pgmat, "pgmat")

        # check xconfig
        check_is_ndarray(xconfig, "xconfig")
        check_ndarray_ndim(xconfig, "xconfig", 2)
        check_ndarray_axis_len(xconfig, "xconfig", 1, self.nparent)

        # check nmating
        if isinstance(nmating, Integral):
            nmating = numpy.repeat(nmating, len(xconfig))
        check_is_ndarray(nmating, "nmating")
        check_ndarray_shape_eq(nmating, "nmating", (len(xconfig),))

        # check nprogeny
        if isinstance(nprogeny, Integral):
            nprogeny = numpy.repeat(nprogeny, len(xconfig))
        check_is_ndarray(nprogeny, "nprogeny")
        check_ndarray_shape_eq(nprogeny, "nprogeny", (len(xconfig),))

        # check miscout
        if miscout is not None:
            check_is_dict(miscout, "miscout")
        
        # check nself
        check_is_Integral(nself, "nself")

        ########################################################################
        ########################## Progeny generation ##########################
        # get female and male selections; repeat by nmating
        fsel = numpy.repeat(xconfig[:,0], nmating * nprogeny)
        msel = numpy.repeat(xconfig[:,1], nmating * nprogeny)

        # get pointers to genotypes and crossover probabilities, respectively
        geno = pgmat.mat
        xoprob = pgmat.vrnt_xoprob

        # create hybrid genotypes
        hgeno = mat_mate(geno, geno, fsel, msel, xoprob, self.rng)

        # generate selection array for all hybrid lines
        asel = numpy.arange(hgeno.shape[1])

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
        nfam = len(xconfig)                     # calculate number of families
        taxa_grp = numpy.repeat(                # construct taxa_grp vector
            numpy.arange(                       # repeat for crosses * progeny
                self.family_counter,            # start family number (inclusive)
                self.family_counter + nfam,     # stop family number (exclusive)
                dtype = 'int64'
            ),
            nmating * nprogeny
        )
        self.family_counter += nfam             # increment counter
        
        ########################################################################
        ########################## Output generation ###########################
        # create new DensePhasedGenotypeMatrix
        progeny = DensePhasedGenotypeMatrix(
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



################################## Utilities ###################################
def check_is_TwoWayCross(v: object, vname: str) -> None:
    """
    Check if object is of type TwoWayCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, TwoWayCross):
        raise TypeError("'%s' must be a TwoWayCross." % vname)
