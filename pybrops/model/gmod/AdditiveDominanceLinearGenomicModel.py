"""
Module defining interfaces and error checking routines for genomic prediction
models that incorporate genomic additive and dominance effects.
"""

from typing import Any
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel

class AdditiveDominanceLinearGenomicModel(AdditiveLinearGenomicModel):
    """
    The AdditiveDominanceLinearGenomicModel class represents an interface for a
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
        \\mathbf{Z} = \\begin{bmatrix} \\mathbf{Z_{misc}} & \\mathbf{Z_{a} & \\mathbf{Z_{d}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{Z_{misc}}` is a matrix of miscellaneous random effect predictors of shape ``(n,p_misc)``
    - :math:`\\mathbf{Z_{a}}` is a matrix of additive genomic marker predictors of shape ``(n,p_a)``
    - :math:`\\mathbf{Z_{d}}` is a matrix of dominance genomic marker predictors of shape ``(n,p_d)``

    .. math::
        \\mathbf{U} = \\begin{bmatrix} \\mathbf{U_{misc}} \\\\ \\mathbf{U_{a}} \\\\ \\mathbf{U_{d}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{U_{misc}}` is a matrix of miscellaneous random effects of shape ``(p_misc,t)``
    - :math:`\\mathbf{U_{a}}` is a matrix of additive genomic marker effects of shape ``(p_a,t)``
    - :math:`\\mathbf{U_{d}}` is a matrix of dominance genomic marker effects of shape ``(p_d,t)``

    Shape definitions:

    - ``n`` is the number of individuals
    - ``q`` is the number of fixed effect predictors (e.g. environments)
    - ``p`` is the number of random effect predictors.
    - ``p_misc`` is the number of miscellaneous random effect predictors.
    - ``p_a`` is the number of additive genomic marker predictors.
    - ``p_d`` is the number of dominance genomic marker predictors.
    - The sum of ``p_misc``, ``p_a``, and ``p_d`` equals ``p``.
    - ``t`` is the number of traits
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class AdditiveDominanceLinearGenomicModel.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(AdditiveDominanceLinearGenomicModel, self).__init__(**kwargs)

    ############################ Object Properties #############################

    @property
    def u_d(self) -> Any:
        """Dominance genomic marker effects."""
        raise NotImplementedError("property is abstract")
    @u_d.setter
    def u_d(self, value: Any) -> None:
        """Set dominance genomic marker effects"""
        raise NotImplementedError("property is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_AdditiveDominanceLinearGenomicModel(v: object, vname: str) -> None:
    """
    Check if object is of type AdditiveDominanceLinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, AdditiveDominanceLinearGenomicModel):
        raise TypeError("variable '{0}' must be a AdditiveDominanceLinearGenomicModel".format(vname))
