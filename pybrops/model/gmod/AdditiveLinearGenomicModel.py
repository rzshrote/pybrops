"""
Module defining interfaces and error checking routines for genomic prediction
models that incorporate genomic additive effects.
"""

from abc import ABCMeta
from abc import abstractmethod
from numbers import Integral
from pybrops.model.gmod.LinearGenomicModel import LinearGenomicModel

class AdditiveLinearGenomicModel(
        LinearGenomicModel,
        metaclass = ABCMeta,
    ):
    """
    The AdditiveLinearGenomicModel class represents an interface for a
    Multivariate Multiple Linear Regression model.

    A Multivariate Multiple Linear Regression model is defined as:

    .. math::
        \\mathbf{Y} = \\mathbf{XB} + \\mathbf{ZU} + \\mathbf{E}

    Where:

    - :math:`\\mathbf{Y}` is a matrix of response variables of shape ``(n,t)``.
    - :math:`\\mathbf{X}` is a matrix of fixed effect predictors of shape ``(n,q)``.
    - :math:`\\mathbf{B}` is a matrix of fixed effect regression coefficients of shape ``(q,t)``.
    - :math:`\\mathbf{Z}` is a matrix of random effect predictors of shape ``(n,p)``.
    - :math:`\\mathbf{U}` is a matrix of random effect regression coefficients of shape ``(p,t)``.
    - :math:`\\mathbf{E}` is a matrix of error terms of shape ``(n,t)``.

    Block matrix modifications to :

    :math:`\\mathbf{Z}` and :math:`\\mathbf{U}` can be decomposed into block
    matrices pertaining to different sets of effects:

    .. math::
        \\mathbf{Z} = \\begin{bmatrix} \\mathbf{Z_{misc}} & \\mathbf{Z_{a}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{Z_{misc}}` is a matrix of miscellaneous random effect predictors of shape ``(n,p_misc)``
    - :math:`\\mathbf{Z_{a}}` is a matrix of additive genomic marker predictors of shape ``(n,p_a)``

    .. math::
        \\mathbf{U} = \\begin{bmatrix} \\mathbf{U_{misc}} \\\\ \\mathbf{U_{a}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{U_{misc}}` is a matrix of miscellaneous random effects of shape ``(p_misc,t)``
    - :math:`\\mathbf{U_{a}}` is a matrix of additive genomic marker effects of shape ``(p_a,t)``

    Shape definitions:

    - ``n`` is the number of individuals
    - ``q`` is the number of fixed effect predictors (e.g. environments)
    - ``p`` is the number of random effect predictors.
    - ``p_misc`` is the number of miscellaneous random effect predictors.
    - ``p_a`` is the number of additive genomic marker predictors.
    - The sum of ``p_misc`` and ``p_a`` equals ``p``.
    - ``t`` is the number of traits
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ####### Additive Linear Genomic Model Parameters #######
    @property
    @abstractmethod
    def nexplan_u_misc(self) -> Integral:
        """Number of miscellaneous random effect explanatory variables required by the model."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def nparam_u_misc(self) -> Integral:
        """Number of miscellaneous random effect parameters."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def u_misc(self) -> object:
        """Miscellaneous random effects."""
        raise NotImplementedError("property is abstract")
    @u_misc.setter
    @abstractmethod
    def u_misc(self, value: object) -> None:
        """Set miscellaneous random effects"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def nexplan_u_a(self) -> Integral:
        """Number of additive genomic marker explanatory variables required by the model."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def nparam_u_a(self) -> Integral:
        """Number of additive genomic marker parameters."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def u_a(self) -> object:
        """Additive genomic marker effects."""
        raise NotImplementedError("property is abstract")
    @u_a.setter
    @abstractmethod
    def u_a(self, value: object) -> None:
        """Set additive genomic marker effects"""
        raise NotImplementedError("property is abstract")



################################## Utilities ###################################
def check_is_AdditiveLinearGenomicModel(v: object, vname: str) -> None:
    """
    Check if object is of type AdditiveLinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, AdditiveLinearGenomicModel):
        raise TypeError("variable '{0}' must be a AdditiveLinearGenomicModel".format(vname))
