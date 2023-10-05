"""
Module defining interfaces and error checking routines for genomic models that
are linear in nature.
"""

from abc import ABCMeta, abstractmethod
from numbers import Integral
from pybrops.core.io.CSVDictInputOutput import CSVDictInputOutput
from pybrops.core.io.PandasDictInputOutput import PandasDictInputOutput
from pybrops.model.gmod.GenomicModel import GenomicModel

class LinearGenomicModel(
        GenomicModel,
        PandasDictInputOutput,
        CSVDictInputOutput,
        metaclass = ABCMeta,
    ):
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

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ########### Linear Genomic Model Parameters ############
    @property
    @abstractmethod
    def nexplan_beta(self) -> Integral:
        """Number of fixed effect explanatory variables required by the model."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def nparam_beta(self) -> Integral:
        """Number of fixed effect parameters."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def beta(self) -> object:
        """Fixed effect regression coefficients."""
        raise NotImplementedError("property is abstract")
    @beta.setter
    @abstractmethod
    def beta(self, value: object) -> None:
        """Set fixed effect regression coefficients"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def nexplan_u(self) -> Integral:
        """Number of random effect explanatory variables required by the model."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def nparam_u(self) -> Integral:
        """Number of random effect parameters."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def u(self) -> object:
        """Random effect regression coefficients."""
        raise NotImplementedError("property is abstract")
    @u.setter
    @abstractmethod
    def u(self, value: object) -> None:
        """Set random effect regression coefficients"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################



################################## Utilities ###################################
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
