"""
Module implementing a dense matrix with axes that are square and associated
error checking routines.
"""

from typing import Any
import numpy

from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_value_numpy import check_ndarray_ndim_gteq
from pybrops.core.mat.DenseMatrix import DenseMatrix
from pybrops.core.mat.SquareMatrix import SquareMatrix

class DenseSquareMatrix(DenseMatrix,SquareMatrix):
    """
    A concrete class for dense square matrices. A "square matrix" is defined as
    a matrix that has the same axis metadata associated with two or more axes.
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

    The purpose of this concrete class is to implement base functionality for:
        1) Square matrix axis metadata.
        2) Determination of square matrix conformity.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DenseSquareMatrix.

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(DenseSquareMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ##################### Matrix Data ######################
    @DenseMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set pointer to raw numpy.ndarray object."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim_gteq(value, "mat", 2)
        self._mat = value

    ############## Square Metadata Properties ##############
    @property
    def nsquare(self) -> int:
        """Number of axes that are square."""
        return len(self.square_axes)
    @nsquare.setter
    def nsquare(self, value: int) -> None:
        """Set the number of axes that are square"""
        error_readonly("nsquare")
    @nsquare.deleter
    def nsquare(self) -> None:
        """Delete the number of axes that are square"""
        error_readonly("nsquare")

    @property
    def square_axes(self) -> tuple:
        """Axis indices for axes that are square."""
        return (0,1)
    @square_axes.setter
    def square_axes(self, value: tuple) -> None:
        """Set axis indices for axes that are square"""
        error_readonly("square_axes")
    @square_axes.deleter
    def square_axes(self) -> None:
        """Delete axis indices for axes that are square"""
        error_readonly("square_axes")

    @property
    def square_axes_len(self) -> tuple:
        """Axis lengths for axes that are square."""
        return tuple(self.mat_shape[ix] for ix in self.square_axes)
    @square_axes_len.setter
    def square_axes_len(self, value: tuple) -> None:
        """Set axis lengths for axes that are square"""
        error_readonly("square_axes_len")
    @square_axes_len.deleter
    def square_axes_len(self) -> None:
        """Delete axis lengths for axes that are square"""
        error_readonly("square_axes_len")

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

    ############################## Object Methods ##############################

    #################### Square Methods ####################
    def is_square(
            self
        ) -> bool:
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
def is_DenseSquareMatrix(v: object) -> bool:
    """
    Determine whether an object is a ``DenseSquareMatrix``.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``DenseSquareMatrix`` object instance.
    """
    return isinstance(v, DenseSquareMatrix)

def check_is_DenseSquareMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``DenseSquareMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseSquareMatrix):
        raise TypeError("'{0}' must be a DenseSquareMatrix".format(vname))
