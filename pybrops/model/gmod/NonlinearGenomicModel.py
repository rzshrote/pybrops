"""
Module defining interfaces and error checking routines for genomic models that
are non-linear in nature.
"""

from abc import ABCMeta
from pybrops.model.gmod.GenomicModel import GenomicModel

class NonlinearGenomicModel(
        GenomicModel,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for non-linear genomic models.

    The purpose for this abstract interface is to provide an interface through
    which non-linear models (e.g. neural networks) may be incorporated into
    PyBrOpS.
    """

    ########################## Special Object Methods ##########################



################################## Utilities ###################################
def check_is_NonlinearGenomicModel(v: object, vname: str) -> None:
    """
    Check if object is of type NonlinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, NonlinearGenomicModel):
        raise TypeError("variable '{0}' must be a NonlinearGenomicModel".format(vname))
