"""
Module defining interfaces and error checking routines for genomic prediction
models that incorporate genomic coancestry effects.
"""

from typing import Any
from pybrops.model.gmod.LinearGenomicModel import LinearGenomicModel

class CoancestryLinearGenomicModel(LinearGenomicModel):
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

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class CoancestryLinearGenomicModel.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(CoancestryLinearGenomicModel, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def u_misc():
        doc = "Miscellaneous random effects."
        def fget(self):
            """Get miscellaneous random effects"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set miscellaneous random effects"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete miscellaneous random effects"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    u_misc = property(**u_misc())

    def u_c():
        doc = "Genomic coancestry effects."
        def fget(self):
            """Get genomic coancestry effects"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set genomic coancestry effects"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete genomic coancestry effects"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    u_c = property(**u_c())



################################################################################
################################## Utilities ###################################
################################################################################
def is_CoancestryLinearGenomicModel(v: Any) -> bool:
    """
    Determine whether an object is a CoancestryLinearGenomicModel.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a CoancestryLinearGenomicModel object instance.
    """
    return isinstance(v, CoancestryLinearGenomicModel)

def check_is_CoancestryLinearGenomicModel(v: Any, vname: str) -> None:
    """
    Check if object is of type CoancestryLinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, CoancestryLinearGenomicModel):
        raise TypeError("variable '{0}' must be a CoancestryLinearGenomicModel".format(vname))
