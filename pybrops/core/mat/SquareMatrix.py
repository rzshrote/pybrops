"""
Module defining interfaces and associated error checking routines for matrices
with axes that are square.
"""

__all__ = [
    "SquareMatrix",
    "check_is_SquareMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from pybrops.core.mat.Matrix import Matrix

class SquareMatrix(
        Matrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for square matrices. A "square matrix" is defined as a
    matrix that has the same axis metadata associated with two or more axes.
    For example::

        This is a square matrix since metadata applies to axes 0 and 1:
               taxa
             +-------+
        taxa | (n,n) |
             +-------+

        This is not a square matrix since unique metadata applies to each axis:
               vrnt
             +-------+
        taxa | (n,p) |  Where: n == p
             +-------+

    The purpose of this abstract class is to provide base functionality for:

        1. Square matrix axis metadata.
        2. Determination of square matrix conformity.
    """

    ########################## Special Object Methods ##########################

    ############## Forward numeric operators ###############
    ### __add__         inherited from ``Matrix``
    ### __sub__         inherited from ``Matrix``
    ### __mul__         inherited from ``Matrix``
    ### __matmul__      inherited from ``Matrix``
    ### __truediv__     inherited from ``Matrix``
    ### __floordiv__    inherited from ``Matrix``
    ### __mod__         inherited from ``Matrix``
    ### __divmod__      inherited from ``Matrix``
    ### __pow__         inherited from ``Matrix``
    ### __lshift__      inherited from ``Matrix``
    ### __rshift__      inherited from ``Matrix``
    ### __and__         inherited from ``Matrix``
    ### __xor__         inherited from ``Matrix``
    ### __or__          inherited from ``Matrix``

    ############# Backwards numeric operators ##############
    ### __radd__        inherited from ``Matrix``
    ### __rsub__        inherited from ``Matrix``
    ### __rmul__        inherited from ``Matrix``
    ### __rmatmul__     inherited from ``Matrix``
    ### __rtruediv__    inherited from ``Matrix``
    ### __rfloordiv__   inherited from ``Matrix``
    ### __rmod__        inherited from ``Matrix``
    ### __rdivmod__     inherited from ``Matrix``
    ### __rlshift__     inherited from ``Matrix``
    ### __rrshift__     inherited from ``Matrix``
    ### __rand__        inherited from ``Matrix``
    ### __rxor__        inherited from ``Matrix``
    ### __ror__         inherited from ``Matrix``

    ############# Augmented numeric operators ##############
    ### __iadd__        inherited from ``Matrix``
    ### __isub__        inherited from ``Matrix``
    ### __imul__        inherited from ``Matrix``
    ### __imatmul__     inherited from ``Matrix``
    ### __itruediv__    inherited from ``Matrix``
    ### __ifloordiv__   inherited from ``Matrix``
    ### __imod__        inherited from ``Matrix``
    ### __ipow__        inherited from ``Matrix``
    ### __ilshift__     inherited from ``Matrix``
    ### __irshift__     inherited from ``Matrix``
    ### __iand__        inherited from ``Matrix``
    ### __ixor__        inherited from ``Matrix``
    ### __ior__         inherited from ``Matrix``

    ################## Logical operators ###################
    ### __lt__          inherited from ``Matrix``
    ### __le__          inherited from ``Matrix``
    ### __eq__          inherited from ``Matrix``
    ### __ne__          inherited from ``Matrix``
    ### __gt__          inherited from ``Matrix``
    ### __ge__          inherited from ``Matrix``

    ################# Container operators ##################
    ### __len__         inherited from ``Matrix``
    ### __getitem__     inherited from ``Matrix``
    ### __setitem__     inherited from ``Matrix``
    ### __delitem__     inherited from ``Matrix``
    ### __iter__        inherited from ``Matrix``

    #################### Matrix copying ####################
    ### __copy__        inherited from ``Matrix``
    ### __deepcopy__    inherited from ``Matrix``

    ########### Miscellaneous special functions ############
    ### __repr__        inherited from ``Matrix``

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    ### mat             inherited from ``Matrix``

    ############## Matrix Metadata Properties ##############
    ### mat_ndim        inherited from ``Matrix``
    ### mat_shape       inherited from ``Matrix``

    ############## Square Metadata Properties ##############
    @property
    @abstractmethod
    def nsquare(self) -> int:
        """Number of axes that are square."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def square_axes(self) -> tuple:
        """Axis indices for axes that are square."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def square_axes_len(self) -> tuple:
        """Axis lengths for axes that are square."""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    ### copy            inherited from ``Matrix``
    ### deepcopy        inherited from ``Matrix``

    ######### Matrix element copy-on-manipulation ##########
    ### adjoin          inherited from ``Matrix``
    ### delete          inherited from ``Matrix``
    ### insert          inherited from ``Matrix``
    ### select          inherited from ``Matrix``

    #################### Square Methods ####################
    @abstractmethod
    def is_square(
            self
        ) -> bool:
        """
        Determine whether the axes lengths for the square axes are identical.

        Returns
        -------
        out : bool
            ``True`` if all square axes are the same length.
            ``False`` if not all square axes are the same length.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat          inherited from ``Matrix``



################################## Utilities ###################################
def check_is_SquareMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``SquareMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, SquareMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,SquareMatrix.__name__,type(v).__name__))
