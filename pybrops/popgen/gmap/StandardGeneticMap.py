"""
Module implementing a standard genetic map format and associated error checking
routines.
"""

from typing import Any

import numpy
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_numpy import check_ndarray_dtype
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.popgen.gmap.GeneticMap import GeneticMap

class StandardGeneticMap(GeneticMap):
    """
    A concrete class for representing a standard genetic map format.

    The purpose of this concrete class is to implement functionality for:
        1) Genetic map representation.
        2) Genetic map metadata.
        3) Genetic map routines.
        4) Genetic map interpolation spline construction.
        5) Genetic map spline interpolation.
        6) Import and export of genetic maps.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for creating a standard genetic map object.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
        vrnt_phypos : numpy.ndarray
        vrnt_genpos : numpy.ndarray
        kwargs : dict
            Additional keyword arguments.
        """
        super(StandardGeneticMap, self).__init__(**kwargs)
        self.vrnt_chrgrp = vrnt_chrgrp
        self.vrnt_phypos = vrnt_phypos
        self.vrnt_genpos = vrnt_genpos
        # TODO: check all lengths equivalent

    def __len__(self):
        """Get the number of markers in the genetic map."""
        return len(self._vrnt_genpos)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def vrnt_chrgrp(self) -> numpy.ndarray:
        """Description for property vrnt_chrgrp."""
        return self._vrnt_chrgrp
    @vrnt_chrgrp.setter
    def vrnt_chrgrp(self, value: numpy.ndarray) -> None:
        """Set data for property vrnt_chrgrp."""
        check_is_ndarray(value, "vrnt_chrgrp")
        check_ndarray_dtype(value, "vrnt_chrgrp", numpy.int64)
        check_ndarray_ndim(value, "vrnt_chrgrp", 1)
        self._vrnt_chrgrp = value
    @vrnt_chrgrp.deleter
    def vrnt_chrgrp(self) -> None:
        """Delete data for property vrnt_chrgrp."""
        del self._vrnt_chrgrp

    @property
    def vrnt_phypos(self) -> numpy.ndarray:
        """Description for property vrnt_phypos."""
        return self._vrnt_phypos
    @vrnt_phypos.setter
    def vrnt_phypos(self, value: numpy.ndarray) -> None:
        """Set data for property vrnt_phypos."""
        check_is_ndarray(value, "vrnt_phypos")
        check_ndarray_dtype(value, "vrnt_phypos", numpy.int64)
        check_ndarray_ndim(value, "vrnt_phypos", 1)
        self._vrnt_phypos = value
    @vrnt_phypos.deleter
    def vrnt_phypos(self) -> None:
        """Delete data for property vrnt_phypos."""
        del self._vrnt_phypos

    @property
    def vrnt_genpos(self) -> numpy.ndarray:
        """Description for property vrnt_genpos."""
        return self._vrnt_genpos
    @vrnt_genpos.setter
    def vrnt_genpos(self, value: numpy.ndarray) -> None:
        """Set data for property vrnt_genpos."""
        check_is_ndarray(value, "vrnt_genpos")
        check_ndarray_dtype(value, "vrnt_genpos", numpy.float64)
        check_ndarray_ndim(value, "vrnt_genpos", 1)
        self._vrnt_genpos = value
    @vrnt_genpos.deleter
    def vrnt_genpos(self) -> None:
        """Delete data for property vrnt_genpos."""
        del self._vrnt_genpos

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_StandardGeneticMap(v: Any) -> bool:
    """
    Determine whether an object is a StandardGeneticMap.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a StandardGeneticMap object instance.
    """
    return isinstance(v, StandardGeneticMap)

def check_is_StandardGeneticMap(v: Any, varname: str) -> None:
    """
    Check if object is of type StandardGeneticMap. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, StandardGeneticMap):
        raise TypeError("'{0}' must be of type StandardGeneticMap.".format(varname))
