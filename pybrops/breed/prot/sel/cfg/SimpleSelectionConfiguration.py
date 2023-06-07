"""
Module containing abstract class definitions for selection configurations
"""

from numbers import Integral
from typing import Union

import numpy
from pybrops.breed.prot.sel.cfg.SelectionConfiguration import SelectionConfiguration
from pybrops.core.error.error_type_numpy import check_is_ndarray, check_ndarray_dtype_is_integer
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_shape_eq
from pybrops.core.error.error_value_python import check_is_gt

from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix


class SimpleSelectionConfiguration(SelectionConfiguration):
    """
    A simple selection configuration class containing the basic necessities for
    a SelectionConfiguration object.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ncross: Integral,
            nparent: Integral,
            nmating: Union[Integral,numpy.ndarray],
            nprogeny: Union[Integral,numpy.ndarray],
            pgmat: PhasedGenotypeMatrix,
            xconfig: numpy.ndarray,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseSelectionConfiguration.
        
        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            A genome matrix containing parental candidates
        xconfig : numpy.ndarray
            A mating configuration matrix of shape ``(ncross,nparent)``.
            This matrix contains indices corresponding to parents in ``pgmat``
            and specifies the manner in which individuals are to be mated.
        ncross : Integral
            Number of cross configurations to consider. Example: ``ncross = 10, nparent = 2``
            specifies 10 two-way crosses.
        nparent : Integral
            The number of parents in a given cross configuration. Example: ``ncross = 10, nparent = 2``
            specifies 10 two-way crosses.
        nmating : Integral, numpy.ndarray
            The number of times an individual cross configuration is executed.
            This becomes important in four-way crosses with heterozygous parents, where
            initial F1 hybrids are unique and can affect the dihybrid composition.
        nprogeny : Integral, numpy.ndarray
            The number of progeny to derive from a mating event.
        """
        # order dependent assignments!
        # set shape parameters first
        self.ncross = ncross
        self.nparent = nparent
        # mating parameters second
        self.nmating = nmating
        self.nprogeny = nprogeny
        # set genotypes and cross configuration third
        self.pgmat = pgmat
        self.xconfig = xconfig

    ############################ Object Properties #############################
    @property
    def ncross(self) -> Integral:
        """Number of cross configurations to consider."""
        return self._ncross
    @ncross.setter
    def ncross(self, value: Integral) -> None:
        """Set number of cross configurations to consider."""
        check_is_Integral(value, "ncross")
        check_is_gt(value, "ncross", 0)
        self._ncross = value
    
    @property
    def nparent(self) -> Integral:
        """Number of parents in a given cross configuration."""
        return self._nparent
    @nparent.setter
    def nparent(self, value: Integral) -> None:
        """Set number of parents in a given cross configuration."""
        check_is_Integral(value, "nparent")
        check_is_gt(value, "nparent", 0)
        self._nparent = value
    
    @property
    def nmating(self) -> numpy.ndarray:
        """Number of times an individual cross configuration is executed."""
        return self._nmating
    @nmating.setter
    def nmating(self, value: Union[Integral,numpy.ndarray]) -> None:
        """Set number of times an individual cross configuration is executed."""
        if isinstance(value, Integral):
            value = numpy.repeat(value, self.ncross)
        elif isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_integer(value, "nmating")
        else:
            raise TypeError("variable 'nmating' must be of type '{0}' but received type '{1}'".format(Integral.__name__,type(value).__name__))
        self._nmating = value
    
    @property
    def nprogeny(self) -> numpy.ndarray:
        """Number of progeny to derive from a mating event."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: Union[Integral,numpy.ndarray]) -> None:
        """Set number of progeny to derive from a mating event."""
        if isinstance(value, Integral):
            value = numpy.repeat(value, self.ncross)
        elif isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_integer(value, "nprogeny")
        else:
            raise TypeError("variable 'nprogeny' must be of type '{0}' but received type '{1}'".format(Integral.__name__,type(value).__name__))
        self._nprogeny = value

    @property
    def pgmat(self) -> PhasedGenotypeMatrix:
        """Genome matrix for the parental candidates."""
        return self._pgmat
    @pgmat.setter
    def pgmat(self, value: PhasedGenotypeMatrix) -> None:
        """Set genome matrix for the parental candidates."""
        check_is_PhasedGenotypeMatrix(value, "pgmat")
        self._pgmat = value
    
    @property
    def xconfig(self) -> numpy.ndarray:
        """xconfig."""
        return self._xconfig
    @xconfig.setter
    def xconfig(self, value: numpy.ndarray) -> None:
        """Set xconfig."""
        check_is_ndarray(value, "xconfig")
        check_ndarray_dtype_is_integer(value, "xconfig")
        check_ndarray_shape_eq(value, "xconfig", (self.ncross,self.nparent))
        self._xconfig = value
    
    ############################## Object Methods ##############################

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################
