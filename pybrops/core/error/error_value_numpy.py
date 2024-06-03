"""
Module containing subroutines to check ``numpy.ndarray`` values.
"""

__all__ = [
    "check_ndarray_has_value",
    "check_ndarray_has_values",
    "check_ndarray_in_interval",
    "check_ndarray_all_gt",
    "check_ndarray_all_gteq",
    "check_ndarray_eq",
    "check_ndarray_is_binary",
    "check_ndarray_ndim",
    "check_ndarray_ndim_gteq",
    "check_ndarray_size",
    "check_ndarray_sum",
    "check_ndarray_mean_is_approx",
    "check_ndarray_std_is_approx",
    "check_ndarray_len_eq",
    "check_ndarray_len_gteq",
    "check_ndarray_shape_eq",
    "check_ndarray_axis_len",
    "check_ndarray_axis_len_eq",
    "check_ndarray_axis_len_lt",
    "check_ndarray_axis_len_lteq",
    "check_ndarray_axis_len_gt",
    "check_ndarray_axis_len_gteq",
    "check_ndarray_is_hypercube",
    "check_ndarray_is_square",
    "check_ndarray_is_square_along_axes",
    "check_ndarray_is_triu",
    "check_ndarray_len_is_multiple_of",
]

from numbers import Real
from typing import Optional, Tuple, Union
import numpy

from pybrops.core.error.error_generic_numpy import generic_check_ndarray_shape

################################################################################
############################### check functions ################################
################################################################################

def check_ndarray_has_value(v: numpy.ndarray, vname: str, value: object) -> None:
    """
    Check if a ``numpy.ndarray`` contains a required value.

    Parameters
    ----------
    v : numpy.ndarray
        Input ``numpy.ndarray``.
    vname : str
        Name of the ``numpy.ndarray`` variable.
    value : object
        A required value that the ``numpy.ndarray`` must contain.
    """
    if value not in v:
        raise ValueError("numpy.ndarray '{0}' must have value '{1}'".format(vname,value))

def check_ndarray_has_values(v: numpy.ndarray, vname: str, *args: Tuple[object,...]) -> None:
    """
    Check if a ``numpy.ndarray`` contains all required values.

    Parameters
    ----------
    v : numpy.ndarray
        Input ``numpy.ndarray``.
    vname : str
        Name of the ``numpy.ndarray`` variable.
    args : tuple
        A tuple of required values that the ``numpy.ndarray`` must contain.
    """
    for value in args:
        if value not in v:
            raise ValueError("numpy.ndarray '{0}' must have value '{1}'".format(vname,value))

def check_ndarray_in_interval(v: numpy.ndarray, vname: str, vmin: Real, vmax: Real) -> None:
    """
    Check whether all values in a ``numpy.ndarray`` are in the interval ``[vmin, vmax]``.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check values.
    vname : str
        Name assigned to the input array.
    vmin : Real
        Lower bound (inclusive).
    vmax : Real
        Upper bound (inclusive).
    """
    if numpy.any(v < vmin) or numpy.any(v > vmax):
        raise ValueError("numpy.ndarray '{0}' must ahve all values in the interval [{1}, {2}]".format(vname, vmin, vmax))

def check_ndarray_all_gt(v: numpy.ndarray, vname: str, vmin: Real) -> None:
    """
    Check if all elements in a ``numpy.ndarray`` are all greater than a provided value.

    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vmin : Real
        Minimium value for elements in the array.
    """
    if numpy.any(v <= vmin):
        raise ValueError("numpy.ndarray '{0}' must must have all values greater than {1}".format(vname,vmin))

def check_ndarray_all_gteq(v: numpy.ndarray, vname: str, vmin: Real) -> None:
    """
    Check if all elements in a ``numpy.ndarray`` are all greater than or equal to a provided value.

    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vmin : Real
        Minimium value for elements in the array.
    """
    if numpy.any(v < vmin):
        raise ValueError("numpy.ndarray '{0}' must must have all values greater than or equal to {1}".format(vname,vmin))

################# check_ndarray_eq #################
def check_ndarray_eq(v: numpy.ndarray, vname: str, w: numpy.ndarray, wname: str) -> None:
    """
    Check if two ``numpy.ndarray`` values are equal to each other.
    
    Parameters
    ----------
    v : numpy.ndarray
        First input array.
    vname : str
        Name assigned to the first input array.
    w : numpy.ndarray
        Second input array.
    wname : str
        Name assigned to the second input array.
    """
    if v.shape != w.shape:
        raise ValueError("numpy.ndarrays '{0}' and '{1}' are not equal: must have same shape".format(vname, wname))
    if numpy.any(v != w):
        raise ValueError("numpy.ndarrays '{0}' and '{1}' are not equal: must have same values at all positions".format(vname,wname))

def check_ndarray_is_binary(v: numpy.ndarray, vname: str) -> None:
    """
    Check if all elements in a ``numpy.ndarray`` are either ``0`` or ``1``.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    """
    if not numpy.all((v == 0) | (v == 1)):
        raise ValueError("numpy.ndarray '{0}' must have values equal to 0 or 1".format(vname))

################ check_ndarray_ndim ################
def check_ndarray_ndim(v: numpy.ndarray, vname: str, vndim: int) -> None:
    """
    Check if a ``numpy.ndarray`` has a specific number of dimensions.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vndim : int
        Expected number of dimensions for the array.
    """
    if v.ndim != vndim:
        raise ValueError("numpy.ndarray '{0}' must have dimension equal to {1}".format(vname, vndim))

def check_ndarray_ndim_gteq(v: numpy.ndarray, vname: str, vndim: int) -> None:
    """
    Check if a ``numpy.ndarray`` has at least a specific number of dimensions.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vndim : int
        Minimum number of dimensions for the array.
    """
    if v.ndim < vndim:
        raise ValueError("numpy.ndarray '{0}' must have dimension greater than or equal to {1}".format(vname, vndim))

################ check_ndarray_size ################
def check_ndarray_size(v: numpy.ndarray, vname: str, vsize: int) -> None:
    """
    Check if a ``numpy.ndarray`` has a specific size.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vsize : int
        Expected size of the array.
    """
    if v.size != vsize:
        raise ValueError("numpy.ndarray '{0}' must have size equal to {1}".format(vname, vsize))

################ check_ndarray_sum #################
def check_ndarray_sum(v: numpy.ndarray, vname: str, vsum: numpy.ndarray, vaxis: int) -> None:
    """
    Check if a ``numpy.ndarray`` has a specific sum along a specific axis.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vsum : numpy.ndarray
        Expected sum of the array.
    vaxis : int
        Axis along which to calculate the sum of the input array.
    """
    if numpy.any(v.sum(vaxis) != vsum):
        raise ValueError("numpy.ndarray '{0}' must have sum equal to {1} along axis {2}".format(vname, vsum, vaxis))

def check_ndarray_mean_is_approx(v: numpy.ndarray, vname: str, vmean: Real, vaxis: int = None, rtol: float = 1e-5, atol: float = 1e-8) -> None:
    """
    Check if a ``numpy.ndarray`` has an approximately equal mean along a specific axis.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vmean : Real
        Expected mean of the array.
    vaxis : int
        Axis along which to calculate the mean of the input array.
    rtol : float
        The relative tolerance parameter.
    atol : float
        The absolute tolerance parameter.
    """
    if not numpy.allclose(v.mean(vaxis), vmean, rtol = rtol, atol = atol):
        raise ValueError("numpy.ndarray '{0}' must have a mean of {1} along axis {2}".format(vname,vmean,vaxis))

def check_ndarray_std_is_approx(v: numpy.ndarray, vname: str, vstd: Real, vaxis: int = None, rtol: float = 1e-5, atol: float = 1e-8) -> None:
    """
    Check if a ``numpy.ndarray`` has an approximately equal standard deviation along a specific axis.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vstd : Real
        Expected standard deviation of the array.
    vaxis : int
        Axis along which to calculate the mean of the input array.
    rtol : float
        The relative tolerance parameter.
    atol : float
        The absolute tolerance parameter.
    """
    if not numpy.allclose(v.std(vaxis), vstd, rtol = rtol, atol = atol):
        raise ValueError("numpy.ndarray '{0}' must have a standard deviation of {1} along axis {2}".format(vname,vstd,vaxis))

################ check_ndarray_len #################
def check_ndarray_len_eq(v: numpy.ndarray, vname: str, vlen: int) -> None:
    """
    Check if a ``numpy.ndarray`` has a specific length.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vlen : int
        Expected length of the input array.
    """
    if len(v) != vlen:
        raise ValueError(
            "numpy.ndarray '{0}' must have length == {1} but received length == {2}".format(vname,vlen,len(v)))

def check_ndarray_len_gteq(v: numpy.ndarray, vname: str, vlen: int) -> None:
    """
    Check if a ``numpy.ndarray`` has at least a specific length.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vlen : int
        Minimum length of the input array.
    """
    if len(v) < vlen:
        raise ValueError("numpy.ndarray '{0}' must have length >= {1} but received length == {2}".format(vname,vlen,len(v)))

############### generic_check_ndarray_shape ################
def check_ndarray_shape_eq(v: numpy.ndarray, vname: str, vshape: tuple) -> None:
    """
    Check if a ``numpy.ndarray`` has a specific shape.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vshape : tuple
        Expected shape of the input array.
    """
    if v.shape != vshape:
        raise ValueError("numpy.ndarray '{0}' must have shape equal to {1}".format(vname, vshape))

def check_ndarray_axis_len(v: numpy.ndarray, vname: str, vaxis: int, vlen: int) -> None:
    """
    Check if a ``numpy.ndarray`` has a specific length along a specific axis.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vaxis : int
        Axis along which to measure the length.
    vlen : int
        Expected length of the input array along the provided axis.
    """
    generic_check_ndarray_shape(v, vname, vlen, vaxis)

def check_ndarray_axis_len_eq(v: numpy.ndarray, vname: str, vaxis: int, vaxislen: int) -> None:
    """
    Check if a ``numpy.ndarray`` has a specific length along a specific axis.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vaxis : int
        Axis along which to measure the length.
    vaxislen : int
        Expected length of the input array along the provided axis.
    """
    if v.shape[vaxis] != vaxislen:
        raise ValueError("numpy.ndarray '{0}' must have axis {1} equal to {2}".format(vname,vaxis,vaxislen))

def check_ndarray_axis_len_lt(v: numpy.ndarray, vname: str, vaxis: int, vaxislen: int) -> None:
    """
    Check if a ``numpy.ndarray`` has at maximum a specific length along a specific axis.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vaxis : int
        Axis along which to measure the length.
    vaxislen : int
        Maximum length of the input array along the provided axis.
    """
    if v.shape[vaxis] >= vaxislen:
        raise ValueError("numpy.ndarray '{0}' must have axis {1} less than {2}".format(vname,vaxis,vaxislen))

def check_ndarray_axis_len_lteq(v: numpy.ndarray, vname: str, vaxis: int, vaxislen: int) -> None:
    """
    Check if a ``numpy.ndarray`` has at maximum a specific length along a specific axis.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vaxis : int
        Axis along which to measure the length.
    vaxislen : int
        Maximum length of the input array along the provided axis.
    """
    if v.shape[vaxis] > vaxislen:
        raise ValueError("numpy.ndarray '{0}' must have axis {1} less than or equal to {2}".format(vname,vaxis,vaxislen))

def check_ndarray_axis_len_gt(v: numpy.ndarray, vname: str, vaxis: int, vaxislen: int) -> None:
    """
    Check if a ``numpy.ndarray`` has at minimum a specific length along a specific axis.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vaxis : int
        Axis along which to measure the length.
    vaxislen : int
        Minimum length of the input array along the provided axis.
    """
    if v.shape[vaxis] <= vaxislen:
        raise ValueError("numpy.ndarray '{0}' must have axis {1} greater than {2}".format(vname,vaxis,vaxislen))

def check_ndarray_axis_len_gteq(v: numpy.ndarray, vname: str, vaxis: int, vaxislen: int) -> None:
    """
    Check if a ``numpy.ndarray`` has at minimum a specific length along a specific axis.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vaxis : int
        Axis along which to measure the length.
    vaxislen : int
        Minimum length of the input array along the provided axis.
    """
    if v.shape[vaxis] < vaxislen:
        raise ValueError("numpy.ndarray '{0}' must have axis {1} greater than or equal to {2}".format(vname,vaxis,vaxislen))

############# generic_check_ndarray_is_square ##############
def check_ndarray_is_hypercube(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` is a hypercube.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    """
    s = v.shape # get shape
    if any(s[0] != e for e in s):
        raise ValueError("numpy.ndarray '{0}' must have equal lengths along all axes".format(vname))

# TODO: replace me with ``check_ndarray_is_square_along_axes``
def check_ndarray_is_square(v: numpy.ndarray, vname: str) -> None:
    """
    Check whether a ``numpy.ndarray`` is square along its last two axes.

    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    """
    if v.ndim < 2:
        raise ValueError("numpy.ndarray '{0}' must have at least two dimensions".format(vname))
    if v.shape[-2] != v.shape[-1]:
        raise ValueError("numpy.ndarray '{0}' must have equal lengths along axes {1} and {2}".format(vname,-2 % v.ndim,-1 % v.ndim))

def check_ndarray_is_square_along_axes(v: numpy.ndarray, vname: str, vaxes: Optional[tuple] = None) -> None:
    """
    Check whether a ``numpy.ndarray`` is square along its axes.

    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    vaxes : tuple, None
        Axes along which to check for square property.
        If ``tuple``, only check along the array dimensions provided in the tuple.
        If ``None``, check along all array dimensions.
    """
    # if the axes are not provided, then use all axes
    if vaxes is None:
        vaxes = tuple(range(v.ndim))
    
    # raise error if not enough axes are provided
    if len(vaxes) < 2:
        raise ValueError("must provide at least 2 axes along which to check, but received {0}: vaxes = {1}".format(len(vaxes),vaxes))
    
    # get the shape of the array
    shape = v.shape

    # get the shape of the square axes
    sshape = tuple(shape[axis] for axis in vaxes)

    # create an iterator for the square axes
    siter = iter(sshape)

    # get the length of the first square axis
    s0 = next(siter)

    # test the remaining axes against the first axis
    if any(s != s0 for s in siter):
        raise ValueError(
            "numpy.ndarray '{0}' must have equal lengths along axes {1}, but received lengths of {2} (full shape: {3})".format(
                vname,vaxes,sshape,shape
            )
        )

def check_ndarray_is_triu(v: numpy.ndarray, vname: str) -> None:
    """
    Check whether a ``numpy.ndarray`` is an upper triangle matrix.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    """
    if numpy.any(v != numpy.triu(v)):
        raise ValueError("numpy.ndarray '{0}' must be an upper triangle matrix".format(vname))

############# check_ndarray_len_is_multiple_of #############
def check_ndarray_len_is_multiple_of(v: numpy.ndarray, vname: str, m: int) -> None:
    """
    Check whether a ``numpy.ndarray`` has a length of a specific multiple.
    
    Parameters
    ----------
    v : numpy.ndarray
        Input array.
    vname : str
        Name assigned to input array.
    m : int
        Length multiple.
    """
    if (len(v) % m) != 0:
        raise ValueError("len({0}) is not a multiple of {1}".format(vname, m))
