"""
Module defining interfaces and error checking routines for genomic models that
are non-linear in nature.
"""

from typing import Any
from pybrops.model.gmod.GenomicModel import GenomicModel

class NonlinearGenomicModel(GenomicModel):
    """
    An abstract class for non-linear genomic models.

    The purpose for this abstract interface is to provide an interface through
    which non-linear models (e.g. neural networks) may be incorporated into
    PyBrOpS.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict) -> None:
        """
        Constructor for NonlinearGenomicModel class.

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(NonlinearGenomicModel, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_NonlinearGenomicModel(v: Any) -> bool:
    """
    Determine whether an object is a NonlinearGenomicModel.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a NonlinearGenomicModel object instance.
    """
    return isinstance(v, NonlinearGenomicModel)

def check_is_NonlinearGenomicModel(v: Any, vname: str) -> None:
    """
    Check if object is of type NonlinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, NonlinearGenomicModel):
        raise TypeError("variable '{0}' must be a NonlinearGenomicModel".format(vname))
