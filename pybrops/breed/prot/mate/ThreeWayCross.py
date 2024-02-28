"""
Module implementing mating protocols for three-way crosses.
"""

from numbers import Integral
from typing import Optional, Union
import numpy
from numpy.random import Generator, RandomState
from pybrops.breed.prot.mate.util import mat_mate
from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState, check_is_ndarray
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_type_python import check_is_Integral, check_is_dict
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len, check_ndarray_ndim, check_ndarray_shape_eq
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import check_DensePhasedGenotypeMatrix_has_vrnt_xoprob, check_is_DensePhasedGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.core.random.prng import global_prng

class ThreeWayCross(MatingProtocol):
    """
    Class implementing mating protocols for three-way crosses.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            progeny_counter: Integral = 0, 
            family_counter: Integral = 0, 
            rng: Union[Generator,RandomState,None] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class ThreeWayCross.

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
        return 3
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
        Mate individuals according to a 3-way mate selection scheme.

        Parameters
        ----------
        pgmat : DensePhasedGenotypeMatrix
            A GenotypeMatrix containing candidate breeding individuals.
        xconfig : numpy.ndarray
            A 1D array of indices of selected individuals of shape ``(k,)``.

            Where:

            - ``k`` is the number of selected individuals.

            Indices are paired as follows:

            - First index is the recurrent parent.
            - Second index is the female parent.
            - Third index is the male parent.

            Example::

                xconfig = [1,5,3,8,2,7]
                recurrent = 1,8
                female = 5,2
                male = 3,7
        nmating : numpy.ndarray
            Number of cross patterns to perform.
        nprogeny : numpy.ndarray
            Number of doubled haploid progeny to generate per cross.
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        s : int, default = 0
            Number of selfing generations post-cross.
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
        # get recurrent, female, and male selections; repeat by ncross
        rsel = numpy.repeat(xconfig[:,0], nmating * nprogeny)  # recurrent parent
        fsel = numpy.repeat(xconfig[:,1], nmating)  # female parent
        msel = numpy.repeat(xconfig[:,2], nmating)  # male parent

        # get pointers to genotypes and crossover probabilities, respectively
        geno = pgmat.mat
        xoprob = pgmat.vrnt_xoprob

        # create F1 genotypes
        f1geno = mat_mate(geno, geno, fsel, msel, xoprob, self.rng)

        # generate selection array for all 2-way hybrid lines
        f1sel = numpy.repeat(
            numpy.arange(f1geno.shape[1]),
            numpy.repeat(nprogeny, nmating)
        )

        # generate three way crosses
        hgeno = mat_mate(geno, f1geno, rsel, f1sel, xoprob, self.rng)

        # get selection array for all 3-way hybrids
        asel = numpy.arange(hgeno.shape[1])

        # self down 3-way hybrids if needed
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
        taxa = numpy.array(["3w"+str(i).zfill(7) for i in riter], dtype = "object")
        self.progeny_counter += progcnt         # increment counter

        # calculate taxa family groupings
        nfam = len(xconfig)                     # calculate number of families
        taxa_grp = numpy.repeat(                # construct taxa_grp
            numpy.repeat(                       # repeat for progeny
                numpy.arange(                   # repeat for crosses
                    self.family_counter,        # start family number (inclusive)
                    self.family_counter + nfam, # stop family number (exclusive)
                    dtype = 'int64'
                ),
                nmating
            ), 
            numpy.repeat(nprogeny, nmating)
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
def check_is_ThreeWayCross(v: object, vname: str) -> None:
    """
    Check if object is of type ThreeWayCross. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, ThreeWayCross):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,ThreeWayCross.__name__,type(v).__name__))
