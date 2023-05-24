"""
Module defining interfaces and error checking routines for genomic models that
are linear in nature.
"""

from typing import Any
from pybrops.model.gmod.GenomicModel import GenomicModel

class LinearGenomicModel(GenomicModel):
    """
    The LinearGenomicModel class represents a Multivariate Multiple Linear
    Regression model.

    A Multivariate Multiple Linear Regression model is defined as:

    .. math::
        Y = X \\Beta + ZU + e

    Where:

    - :math:`Y` is a matrix of response variables of shape ``(n,t)``.
    - :math:`X` is a matrix of fixed effect predictors of shape ``(n,q)``.
    - :math:`\\Beta` is a matrix of fixed effect regression coefficients of
      shape ``(q,t)``.
    - :math:`Z` is a matrix of random effect predictors of shape ``(n,p)``.
    - :math:`U` is a matrix of random effect regression coefficients of shape
      ``(p,t)``.
    - :math:`e` is a matrix of error terms of shape ``(n,t)``.

    Shape definitions:

    - ``n`` is the number of individuals
    - ``q`` is the number of fixed effect predictors (e.g. environments)
    - ``p`` is the number of random effect predictors (e.g. genomic markers)
    - ``t`` is the number of traits
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for LinearGenomicModel class.

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(LinearGenomicModel, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Linear Genomic Model Data ###############
    @property
    def beta(self) -> Any:
        """Fixed effect regression coefficients."""
        raise NotImplementedError("property is abstract")
    @beta.setter
    def beta(self, value: Any) -> None:
        """Set fixed effect regression coefficients"""
        raise NotImplementedError("property is abstract")
    @beta.deleter
    def beta(self) -> None:
        """Delete fixed effect regression coefficients"""
        raise NotImplementedError("property is abstract")

    @property
    def u(self) -> Any:
        """Random effect regression coefficients."""
        raise NotImplementedError("property is abstract")
    @u.setter
    def u(self, value: Any) -> None:
        """Set random effect regression coefficients"""
        raise NotImplementedError("property is abstract")
    @u.deleter
    def u(self) -> None:
        """Delete random effect regression coefficients"""
        raise NotImplementedError("property is abstract")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_LinearGenomicModel(v: object) -> bool:
    """
    Determine whether an object is a LinearGenomicModel.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a LinearGenomicModel object instance.
    """
    return isinstance(v, LinearGenomicModel)

def check_is_LinearGenomicModel(v: object, vname: str) -> None:
    """
    Check if object is of type LinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, LinearGenomicModel):
        raise TypeError("variable '{0}' must be a LinearGenomicModel".format(vname))
