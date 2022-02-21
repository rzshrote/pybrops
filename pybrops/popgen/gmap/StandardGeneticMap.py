"""
Module implementing a standard genetic map format and associated error checking
routines.
"""

from pybrops.core.error import check_is_ndarray
from pybrops.core.error import check_ndarray_dtype
from pybrops.core.error import check_ndarray_ndim
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
    def __init__(self, vrnt_chrgrp, vrnt_phypos, vrnt_genpos, **kwargs):
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
    def vrnt_chrgrp():
        doc = "The vrnt_chrgrp property."
        def fget(self):
            return self._vrnt_chrgrp
        def fset(self, value):
            check_is_ndarray(value, "vrnt_chrgrp")
            check_ndarray_dtype(value, "vrnt_chrgrp", numpy.int64)
            check_ndarray_ndim(value, "vrnt_chrgrp", 1)
            self._vrnt_chrgrp = value
        def fdel(self):
            del self._vrnt_chrgrp
        return locals()
    vrnt_chrgrp = property(**vrnt_chrgrp())

    def vrnt_phypos():
        doc = "The vrnt_phypos property."
        def fget(self):
            return self._vrnt_phypos
        def fset(self, value):
            check_is_ndarray(value, "vrnt_phypos")
            check_ndarray_dtype(value, "vrnt_phypos", numpy.int64)
            check_ndarray_ndim(value, "vrnt_phypos", 1)
            self._vrnt_phypos = value
        def fdel(self):
            del self._vrnt_phypos
        return locals()
    vrnt_phypos = property(**vrnt_phypos())

    def vrnt_genpos():
        doc = "The vrnt_genpos property."
        def fget(self):
            return self._vrnt_genpos
        def fset(self, value):
            check_is_ndarray(value, "vrnt_genpos")
            check_ndarray_dtype(value, "vrnt_genpos", numpy.float64)
            check_ndarray_ndim(value, "vrnt_genpos", 1)
            self._vrnt_genpos = value
        def fdel(self):
            del self._vrnt_genpos
        return locals()
    vrnt_genpos = property(**vrnt_genpos())

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_StandardGeneticMap(v):
    """
    Determine whether an object is a StandardGeneticMap.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a StandardGeneticMap object instance.
    """
    return isinstance(v, StandardGeneticMap)

def check_is_StandardGeneticMap(v, varname):
    """
    Check if object is of type StandardGeneticMap. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_StandardGeneticMap(v):
        raise TypeError("'{0}' must be of type StandardGeneticMap.".format(varname))

def cond_check_is_StandardGeneticMap(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type StandardGeneticMap. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a StandardGeneticMap.
    """
    if cond(v):
        check_is_StandardGeneticMap(v, varname)
