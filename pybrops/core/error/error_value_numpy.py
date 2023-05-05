from typing import Any
import numpy

from . import generic_check_ndarray_shape

################################################################################
############################### check functions ################################
################################################################################

def check_ndarray_in_interval(v: numpy.ndarray, vname: str, vmin: float, vmax: float) -> None:
    if numpy.any(v < vmin) or numpy.any(v > vmax):
        raise ValueError("variable '{0}' is not in interval [{1}, {2}]".format(vname, vmin, vmax))

################# check_ndarray_eq #################
def check_ndarray_eq(v: numpy.ndarray, vname: str, w: numpy.ndarray, wname: str) -> None:
    if v.shape != w.shape:
        raise ValueError("arrays '{0}' and '{1}' are not equal: must have same shape".format(vname, wname))
    if numpy.any(v != w):
        raise ValueError("arrays '{0}' and '{1}' are not equal: must have same values".format(vname,wname))

def check_ndarray_align(v: numpy.ndarray, vname: str, w: numpy.ndarray, wname: str) -> None:
    if v.shape != w.shape:
        raise ValueError("arrays '{0}' and '{1}' are not aligned: must have same shape".format(vname, wname))
    if numpy.any(v != w):
        raise ValueError("arrays '{0}' and '{1}' are not aligned: must have same values at same position".format(vname,wname))

def check_ndarray_is_binary(v: numpy.ndarray, vname: str) -> None:
    if not numpy.all((v == 0) | (v == 1)):
        raise ValueError("variable '{0}' must have values equal to 0 or 1".format(vname))

def check_ndarray_is_positive(v: numpy.ndarray, vname: str) -> None:
    if numpy.any(v < 0):
        raise ValueError("variable '{0}' must have all positive values".format(vname))

################ check_ndarray_ndim ################
def check_ndarray_ndim(v: numpy.ndarray, vname: str, vndim: int) -> None:
    if v.ndim != vndim:
        raise ValueError("variable '{0}' must have dimension equal to {1}".format(vname, vndim))

def check_ndarray_ndim_gteq(v: numpy.ndarray, vname: str, vndim: int) -> None:
    if v.ndim < vndim:
        raise ValueError("variable '{0}' must have dimension greater than or equal to {1}".format(vname, vndim))

################ check_ndarray_size ################
def check_ndarray_size(v: numpy.ndarray, vname: str, vsize: int):
    if v.size != vsize:
        raise ValueError("variable '{0}' must have size equal to {1}".format(vname, vsize))

################ check_ndarray_sum #################
def check_ndarray_sum(v: numpy.ndarray, vname: str, vsum: numpy.ndarray, vaxis: int):
    if numpy.any(v.sum(vaxis) != vsum):
        raise ValueError("variable '{0}' must have sum equal to {1} along axis {2}".format(vname, vsum, vaxis))

def check_ndarray_mean_is_approx(v: numpy.ndarray, vname: str, vmean: float, vaxis: int = None, rtol: float = 1e-5, atol: float = 1e-8):
    if not numpy.allclose(v.mean(vaxis), vmean, rtol = rtol, atol = atol):
        raise ValueError("'{0}' must have a mean of {1} along axis {2}".format(vname,vmean,vaxis))

def check_ndarray_std_is_approx(v, vname, vstd, vaxis = None, rtol = 1e-5, atol = 1e-8):
    if not numpy.allclose(v.std(vaxis), vstd, rtol = rtol, atol = atol):
        raise ValueError("'{0}' must have a standard deviation of {1} along axis {2}".format(vname,vstd,vaxis))

################ check_ndarray_len #################
def check_ndarray_len_eq(v: numpy.ndarray, vname: str, vlen: int) -> None:
    if len(v) != vlen:
        raise ValueError("variable '{0}' must have length equal to {1}".format(vname, vlen))

def check_ndarray_len_gteq(v: numpy.ndarray, vname: str, vlen: int) -> None:
    if len(v) < vlen:
        raise ValueError("variable '{0}' must have length greater than or equal to {1}".format(vname, vlen))

############### generic_check_ndarray_shape ################
def check_ndarray_shape_eq(v: numpy.ndarray, vname: str, vshape: tuple):
    if v.shape != vshape:
        raise ValueError("variable '{0}' must have shape equal to {1}".format(vname, vshape))

def check_ndarray_axis_len(v, vname, vaxis, vlen):
    generic_check_ndarray_shape(v, vname, vlen, vaxis)

def check_ndarray_axis_len_eq(v: numpy.ndarray, vname: str, vaxis: int, vaxislen: int):
    if v.shape[vaxis] != vaxislen:
        raise ValueError("variable '{0}' must have axis {1} equal to {2}".format(vname,vaxis,vaxislen))

def check_ndarray_axis_len_lt(v: numpy.ndarray, vname: str, vaxis: int, vaxislen: int):
    if v.shape[vaxis] >= vaxislen:
        raise ValueError("variable '{0}' must have axis {1} less than {2}".format(vname,vaxis,vaxislen))

def check_ndarray_axis_len_lteq(v: numpy.ndarray, vname: str, vaxis: int, vaxislen: int):
    if v.shape[vaxis] > vaxislen:
        raise ValueError("variable '{0}' must have axis {1} less than or equal to {2}".format(vname,vaxis,vaxislen))

def check_ndarray_axis_len_gt(v: numpy.ndarray, vname: str, vaxis: int, vaxislen: int):
    if v.shape[vaxis] <= vaxislen:
        raise ValueError("variable '{0}' must have axis {1} greater than {2}".format(vname,vaxis,vaxislen))

def check_ndarray_axis_len_gteq(v: numpy.ndarray, vname: str, vaxis: int, vaxislen: int):
    if v.shape[vaxis] < vaxislen:
        raise ValueError("variable '{0}' must have axis {1} greater than or equal to {2}".format(vname,vaxis,vaxislen))

############# generic_check_ndarray_is_square ##############
def check_ndarray_is_square(v: numpy.ndarray, vname: str) -> None:
    s = v.shape # get shape
    if any(s[0] != e for e in s):
        raise ValueError("variable '{0}' must have equal lengths along all axes".format(vname))

def check_ndarray_is_triu(v: numpy.ndarray, vname: str) -> None:
    if numpy.any(v != numpy.triu(v)):
        raise ValueError("variable '{0}' must be an upper triangle matrix".format(vname))

############# check_ndarray_len_is_multiple_of #############
def check_ndarray_len_is_multiple_of(v, vname, m):
    if (len(v) % m) != 0:
        raise ValueError("len({0}) is not a multiple of {1}".format(vname, m))
