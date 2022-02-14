"""
Module for providing abstract interfaces and functions pertaining to square
matrices.
"""

from pybrops.core.mat.Matrix import Matrix

class SquareMatrix(Matrix):
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
        1) Square matrix axis metadata.
        2) Determination of square matrix conformity.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the SquareMatrix abstract class.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments. Used for cooperative inheritance.
            Dictionary passing unused arguments to the parent class constructor.
        """
        super(SquareMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Square Metadata Properties ##############
    def nsquare():
        doc = "Number of axes that are square"
        def fget(self):
            """Get the number of axes that are square"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the number of axes that are square"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the number of axes that are square"""
            raise NotImplementedError("method is abstract")
        return locals()
    nsquare = property(**nsquare())

    def square_axes():
        doc = "Axis indices for axes that are square"
        def fget(self):
            """Get axis indices for axes that are square"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set axis indices for axes that are square"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete axis indices for axes that are square"""
            raise NotImplementedError("method is abstract")
        return locals()
    square_axes = property(**square_axes())

    def square_axes_len():
        doc = "Axis lengths for axes that are square"
        def fget(self):
            """Get axis lengths for axes that are square"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set axis lengths for axes that are square"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete axis lengths for axes that are square"""
            raise NotImplementedError("method is abstract")
        return locals()
    square_axes_len = property(**square_axes_len())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    #################### Square Methods ####################
    def is_square(self):
        """
        Determine whether the axis lengths for the square axes are identical.

        Returns
        -------
        out : bool
            ``True`` if all square axes are the same length.
            ``False`` if not all square axes are the same length.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_SquareMatrix(obj):
    """
    Determine whether an object is a ``SquareMatrix``.

    Parameters
    ----------
    obj : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``SquareMatrix`` object instance.
    """
    return isinstance(obj, SquareMatrix)

def check_is_SquareMatrix(obj, objname):
    """
    Check if object is of type ``SquareMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(obj, SquareMatrix):
        raise TypeError("'{0}' must be a SquareMatrix".format(objname))

def cond_check_is_SquareMatrix(obj, objname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type ``SquareMatrix``. Otherwise raise
    ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    cond : function
        A function returning ``True`` or ``False`` for whether to test if ``obj``
        is a ``SquareMatrix``.
    """
    if cond(obj) and not isinstance(obj, SquareMatrix):
        raise TypeError("'{0}' must be a SquareMatrix".format(objname))
