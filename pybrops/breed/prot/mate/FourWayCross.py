"""
Module implementing mating protocols for four-way crosses.
"""

from numbers import Integral
from typing import Any, Optional, Union
import numpy

from pybrops.breed.prot.mate.util import mat_mate
from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_type_python import check_is_int
from pybrops.popgen.gmat import PhasedGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import check_is_DensePhasedGenotypeMatrix
from pybrops.core.random.prng import global_prng

class FourWayCross(MatingProtocol):
    """
    Class implementing mating protocols for four-way crosses.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            progeny_counter: int = 0, 
            family_counter: int = 0, 
            rng: Union[numpy.random.Generator,numpy.random.RandomState,None] = global_prng, 
            **kwargs: dict
        ) -> None:
        super(FourWayCross, self).__init__(**kwargs)

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
        return 4
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
            xconfig: numpy.ndarray, 
            ncross: Union[int,numpy.ndarray], 
            nprogeny: Union[int,numpy.ndarray], 
            miscout: Optional[dict] = None, 
            nself: int = 0, 
            **kwargs: dict
        ) -> PhasedGenotypeMatrix:
        """
        Mate individuals according to a 4-way mate selection scheme.

        Parameters
        ----------
        pgmat : DensePhasedGenotypeMatrix
            A GenotypeMatrix containing candidate breeding individuals.
        xconfig : numpy.ndarray
            A 1D array of indices of selected individuals of shape ``(k,)``.

            Where:

            - ``k`` is the number of selected individuals.

            Indices are paired as follows:

            - First index is the female parent 2.
            - Second index is the male parent 2.
            - Third index is the female parent 1.
            - Fourth index is the male parent 1.

            Example::

                xconfig = [1,5,3,8]
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
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        s : int, default = 0
            Number of selfing generations post-cross.
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

        # get female2, male2, female1, and male1 selections; repeat by ncross
        f2sel = numpy.repeat(xconfig[0::4], ncross)
        m2sel = numpy.repeat(xconfig[1::4], ncross)
        f1sel = numpy.repeat(xconfig[2::4], ncross)
        m1sel = numpy.repeat(xconfig[3::4], ncross)

        # get pointers to genotypes and crossover probabilities, respectively
        geno = pgmat.mat
        xoprob = pgmat.vrnt_xoprob

        # create F1 genotypes
        abgeno = mat_mate(geno, geno, f1sel, m1sel, xoprob, self.rng)
        cdgeno = mat_mate(geno, geno, f2sel, m2sel, xoprob, self.rng)

        # generate selection array for all hybrid lines
        asel = numpy.repeat(numpy.arange(abgeno.shape[1]), nprogeny)

        # generate dihybrid cross
        hgeno = mat_mate(abgeno, cdgeno, asel, asel, xoprob, self.rng)

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
        taxa = numpy.array(["3w"+str(i).zfill(7) for i in riter], dtype = "object")
        self.progeny_counter += progcnt         # increment counter

        # calculate taxa family groupings
        nfam = len(xconfig) // 4                    # calculate number of families
        taxa_grp = numpy.repeat(                # construct taxa_grp
            numpy.repeat(                       # repeat for progeny
                numpy.arange(                   # repeat for crosses
                    self.family_counter,        # start family number (inclusive)
                    self.family_counter + nfam, # stop family number (exclusive)
                    dtype = 'int64'
                ),
                ncross
            ), 
            nprogeny
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
def is_FourWayCross(v: object) -> bool:
    return isinstance(v, FourWayCross)

def check_is_FourWayCross(v: object, vname: str) -> None:
    if not isinstance(v, FourWayCross):
        raise TypeError("'%s' must be a FourWayCross." % vname)
