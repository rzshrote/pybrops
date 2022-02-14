import numpy

from pybrops.core.mat.SquareMatrix import SquareMatrix
from pybrops.core.mat.DenseMatrix import DenseMatrix

from pybrops.core.error import check_is_ndarray
from pybrops.core.error import check_ndarray_at_least_2d
from pybrops.core.error import error_readonly

class DenseSquareMatrix(DenseMatrix,SquareMatrix):
    """docstring for DenseSquareMatrix."""

    def __init__(self, mat, **kwargs):
        super(DenseSquareMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ##################### Matrix Data ######################
    def mat():
        doc = "The mat property."
        def fget(self):
            return self._mat
        def fset(self, value):
            check_is_ndarray(value, "mat")
            check_ndarray_at_least_2d(value, "mat")
            self._mat = value
        def fdel(self):
            del self._mat
        return locals()
    mat = property(**mat())

    ############## Square Metadata Properties ##############
    def nsquare():
        doc = "Number of axes that are square"
        def fget(self):
            """Get the number of axes that are square"""
            return len(self.square_axes)
        def fset(self, value):
            """Set the number of axes that are square"""
            error_readonly("nsquare")
        def fdel(self):
            """Delete the number of axes that are square"""
            error_readonly("nsquare")
        return locals()
    nsquare = property(**nsquare())

    def square_axes():
        doc = "Axis indices for axes that are square"
        def fget(self):
            """Get axis indices for axes that are square"""
            return (0, 1)
        def fset(self, value):
            """Set axis indices for axes that are square"""
            error_readonly("square_axes")
        def fdel(self):
            """Delete axis indices for axes that are square"""
            error_readonly("square_axes")
        return locals()
    square_axes = property(**square_axes())

    def square_axes_len():
        doc = "Axis lengths for axes that are square"
        def fget(self):
            """Get axis lengths for axes that are square"""
            return tuple(self.mat_shape[ix] for ix in self.square_axes)
        def fset(self, value):
            """Set axis lengths for axes that are square"""
            error_readonly("square_axes_len")
        def fdel(self):
            """Delete axis lengths for axes that are square"""
            error_readonly("square_axes_len")
        return locals()
    square_axes_len = property(**square_axes_len())

    ################### Fill data lookup ###################
    # map dtypes to fill values
    # fill values are the most extreme values from zero if integer, or NaN if floating
    _fill_value = {
        # strings
        "int8": numpy.int8(numpy.iinfo("int8").min),        # -128
        "int16": numpy.int16(numpy.iinfo("int16").min),     # -32768
        "int32": numpy.int32(numpy.iinfo("int32").min),     # -2147483648
        "int64": numpy.int64(numpy.iinfo("int64").min),     # -9223372036854775808
        "int": numpy.int0(numpy.iinfo("int").min),          # system dependent
        "int0": numpy.int0(numpy.iinfo("int0").min),        # system dependent

        "uint8": numpy.uint8(numpy.iinfo("uint8").max),     # 255
        "uint16": numpy.uint16(numpy.iinfo("uint16").max),  # 65535
        "uint32": numpy.uint32(numpy.iinfo("uint32").max),  # 4294967295
        "uint64": numpy.uint64(numpy.iinfo("uint64").max),  # 18446744073709551615
        "uint": numpy.uint0(numpy.iinfo("uint").max),       # system dependent
        "uint0": numpy.uint0(numpy.iinfo("uint0").max),     # system dependent

        "float16": numpy.float16(numpy.nan),                # NaN
        "float32": numpy.float32(numpy.nan),                # NaN
        "float64": numpy.float64(numpy.nan),                # NaN
        "float128": numpy.float128(numpy.nan),              # NaN
        "float": numpy.float_(numpy.nan),                   # NaN

        # actual dtypes
        numpy.dtype("int8"): numpy.int8(numpy.iinfo("int8").min),       # -128
        numpy.dtype("int16"): numpy.int16(numpy.iinfo("int16").min),    # -32768
        numpy.dtype("int32"): numpy.int32(numpy.iinfo("int32").min),    # -2147483648
        numpy.dtype("int64"): numpy.int64(numpy.iinfo("int64").min),    # -9223372036854775808
        numpy.dtype("int"): numpy.int0(numpy.iinfo("int").min),         # system dependent
        numpy.dtype("int0"): numpy.int0(numpy.iinfo("int0").min),       # system dependent

        numpy.dtype("uint8"): numpy.uint8(numpy.iinfo("uint8").max),    # 255
        numpy.dtype("uint16"): numpy.uint16(numpy.iinfo("uint16").max), # 65535
        numpy.dtype("uint32"): numpy.uint32(numpy.iinfo("uint32").max), # 4294967295
        numpy.dtype("uint64"): numpy.uint64(numpy.iinfo("uint64").max), # 18446744073709551615
        numpy.dtype("uint"): numpy.uint0(numpy.iinfo("uint").max),      # system dependent
        numpy.dtype("uint0"): numpy.uint0(numpy.iinfo("uint0").max),    # system dependent

        numpy.dtype("float16"): numpy.float16(numpy.nan),               # NaN
        numpy.dtype("float32"): numpy.float32(numpy.nan),               # NaN
        numpy.dtype("float64"): numpy.float64(numpy.nan),               # NaN
        numpy.dtype("float128"): numpy.float128(numpy.nan),             # NaN
        numpy.dtype("float"): numpy.float_(numpy.nan),                  # NaN
    }

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
        siter = iter(self.square_axes_len)  # get iterator of axis lengths
        try:                                # try next
            e0 = next(siter)                # get the first element of siter
        except StopIteration:               # catch StopIteration exception
            return False                    # something is horribly wrong; no axes!
        out = all(e0 == e for e in siter)   # determine if all are equal to e0
        return out



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseSquareMatrix(obj):
    """
    Determine whether an object is a ``DenseSquareMatrix``.

    Parameters
    ----------
    obj : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``DenseSquareMatrix`` object instance.
    """
    return isinstance(obj, DenseSquareMatrix)

def check_is_DenseSquareMatrix(obj, objname):
    """
    Check if object is of type ``DenseSquareMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(obj, DenseSquareMatrix):
        raise TypeError("'{0}' must be a DenseSquareMatrix".format(objname))

def cond_check_is_DenseSquareMatrix(obj, objname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type ``DenseSquareMatrix``. Otherwise raise
    ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    cond : function
        A function returning ``True`` or ``False`` for whether to test if ``obj``
        is a ``DenseSquareMatrix``.
    """
    if cond(obj) and not isinstance(obj, DenseSquareMatrix):
        raise TypeError("'{0}' must be a DenseSquareMatrix".format(objname))
