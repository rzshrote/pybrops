"""
Module defining interfaces and error checking routines for genomic prediction
models that incorporate genomic coancestry effects.
"""

from abc import ABCMeta
from abc import abstractmethod
from pybrops.model.gmod.LinearGenomicModel import LinearGenomicModel

class CoancestryLinearGenomicModel(
        LinearGenomicModel,
        metaclass = ABCMeta,
    ):
    """
    The CoancestryLinearGenomicModel class represents an interface for a
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
        \\mathbf{Z} = \\begin{bmatrix} \\mathbf{Z_{misc}} & \\mathbf{Z_{c}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{Z_{misc}}` is a matrix of miscellaneous random effect predictors of shape ``(n,p_misc)``
    - :math:`\\mathbf{Z_{c}}` is a matrix of coancestry predictors of shape ``(n,p_c)``

    .. math::
        \\mathbf{U} = \\begin{bmatrix} \\mathbf{U_{misc}} \\\\ \\mathbf{U_{c}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{U_{misc}}` is a matrix of miscellaneous random effects of shape ``(p_misc,t)``
    - :math:`\\mathbf{U_{c}}` is a matrix of coancestry effects of shape ``(p_c,t)``

    Shape definitions:

    - ``n`` is the number of individuals
    - ``q`` is the number of fixed effect predictors (e.g. environments)
    - ``p`` is the number of random effect predictors.
    - ``p_misc`` is the number of miscellaneous random effect predictors.
    - ``p_c`` is the number of coancestry predictors.
    - The sum of ``p_misc`` and ``p_c`` equals ``p``.
    - ``t`` is the number of traits
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################
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
    def u_c(self) -> object:
        """Genomic coancestry effects."""
        raise NotImplementedError("property is abstract")
    @u_c.setter
    @abstractmethod
    def u_c(self, value: object) -> None:
        """Set genomic coancestry effects"""
        raise NotImplementedError("property is abstract")



################################## Utilities ###################################
def check_is_CoancestryLinearGenomicModel(v: object, vname: str) -> None:
    """
    Check if object is of type CoancestryLinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, CoancestryLinearGenomicModel):
        raise TypeError("variable '{0}' must be a CoancestryLinearGenomicModel".format(vname))
