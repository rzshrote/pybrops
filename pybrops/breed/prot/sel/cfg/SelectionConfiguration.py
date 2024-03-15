"""
Module containing abstract class definitions for selection configurations
"""

from abc import ABCMeta, abstractmethod
from numbers import Integral
from typing import Union

import numpy
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.core.error.error_type_numpy import check_is_ndarray, check_ndarray_dtype_is_integer
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_all_gt, check_ndarray_len_eq, check_ndarray_ndim, check_ndarray_shape_eq
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix


class SelectionConfiguration(
        metaclass = ABCMeta,
    ):
    """
    docstring for SelectionConfiguration.
    """

    ########################## Special Object Methods ##########################
    # do not define __init__ since this is a semi-abstract class

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
            check_is_gt(value, "nmating", 0)
            value = numpy.repeat(value, self.ncross)
        elif isinstance(value, numpy.ndarray):
            pass
        else:
            raise TypeError("variable 'nmating' must be of type '{0}' or '{1}' but received type '{2}'".format(Integral.__name__,numpy.ndarray.__name__,type(value).__name__))
        check_ndarray_ndim(value, "nmating", 1)
        check_ndarray_len_eq(value, "nmating", self.ncross)
        check_ndarray_dtype_is_integer(value, "nmating")
        check_ndarray_all_gt(value, "nmating", 0)
        self._nmating = value
    
    @property
    def nprogeny(self) -> numpy.ndarray:
        """Number of progeny to derive from a mating event."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: Union[Integral,numpy.ndarray]) -> None:
        """Set number of progeny to derive from a mating event."""
        if isinstance(value, Integral):
            check_is_gt(value, "nprogeny", 0)
            value = numpy.repeat(value, self.ncross)
        elif isinstance(value, numpy.ndarray):
            pass
        else:
            raise TypeError("variable 'nprogeny' must be of type '{0}' or '{1}' but received type '{2}'".format(Integral.__name__,numpy.ndarray.__name__,type(value).__name__))
        check_ndarray_ndim(value, "nprogeny", 1)
        check_ndarray_len_eq(value, "nprogeny", self.ncross)
        check_ndarray_dtype_is_integer(value, "nprogeny")
        check_ndarray_all_gt(value, "nprogeny", 0)
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



################################## Utilities ###################################
def check_is_SelectionConfiguration(v: object, vname: str) -> None:
    """
    Check if object is of type SelectionConfiguration, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SelectionConfiguration):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,SelectionConfiguration.__name__,type(v).__name__))
