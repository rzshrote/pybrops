from typing import Any
import numpy

from . import generic_check_ndarray_eq
from . import generic_check_ndarray_sum
from . import generic_check_ndarray_ndim
from . import generic_check_ndarray_size
from . import generic_check_ndarray_shape
from . import generic_check_ndarray_is_square
from . import generic_check_ndarray_ndim_gteq

################################################################################
############################### check functions ################################
################################################################################

def check_ndarray_in_interval(v: numpy.ndarray, vname: str, vmin: float, vmax: float):
    if numpy.any(v < vmin) or numpy.any(v > vmax):
        raise ValueError(
            "variable '{0}' is not in interval [{1}, {2}]".format(vname, vmin, vmax)
        )

################# generic_check_ndarray_eq #################
def check_ndarray_eq(v: numpy.ndarray, vname: str, w: numpy.ndarray, wname: str):
    if not numpy.all(v == w):
        raise ValueError("variable '{0}' must have values equal to {1}".format(vname, wname))

def check_ndarray_is_binary(v: numpy.ndarray, vname: str):
    if not numpy.all((v == 0) | (v == 1)):
        raise ValueError("variable '{0}' must have values equal to 0 or 1".format(vname))

def check_ndarray_is_positive(v: numpy.ndarray, vname: str):
    if numpy.any(v < 0):
        raise ValueError("variable '{0}' must have all positive values".format(vname))

################ generic_check_ndarray_ndim ################
def check_ndarray_ndim(v: numpy.ndarray, vname: str, vndim: int):
    if v.ndim != vndim:
        raise ValueError("variable '{0}' must have dimension equal to {1}".format(vname, vndim))

def check_ndarray_ndim_gteq(v, vname, vndim):
    generic_check_ndarray_ndim_gteq(v, vname, vndim)

def check_ndarray_is_1d(v: Any, vname: str) -> None:
    generic_check_ndarray_ndim(v, vname, 1)

def check_ndarray_at_least_1d(v: Any, vname: str) -> None:
    generic_check_ndarray_ndim_gteq(v, vname, 1)

def check_ndarray_is_2d(v: Any, vname: str) -> None:
    generic_check_ndarray_ndim(v, vname, 2)

def check_ndarray_at_least_2d(v: Any, vname: str) -> None:
    generic_check_ndarray_ndim_gteq(v, vname, 2)

def check_ndarray_is_3d(v: Any, vname: str) -> None:
    generic_check_ndarray_ndim(v, vname, 3)

def check_ndarray_at_least_3d(v: Any, vname: str) -> None:
    generic_check_ndarray_ndim_gteq(v, vname, 3)

################ generic_check_ndarray_size ################
def check_ndarray_size(v, vname, vsize):
    generic_check_ndarray_size(v, vname, vsize)

################ generic_check_ndarray_sum #################
def check_ndarray_sum(v, vname, vsum, vaxis):
    generic_check_ndarray_sum(v, vname, vsum, vaxis)

def check_ndarray_mean_is_approx(v: numpy.ndarray, vname: str, vmean: float, vaxis: int = None, rtol: float = 1e-5, atol: float = 1e-8):
    if not numpy.allclose(v.mean(vaxis), vmean, rtol = rtol, atol = atol):
        raise ValueError("'{0}' must have a mean of {1} along axis {2}".format(vname,vmean,vaxis))

def check_ndarray_std_is_approx(v, vname, vstd, vaxis = None, rtol = 1e-5, atol = 1e-8):
    if not numpy.allclose(v.std(vaxis), vstd, rtol = rtol, atol = atol):
        raise ValueError("'{0}' must have a standard deviation of {1} along axis {2}".format(vname,vstd,vaxis))

############### generic_check_ndarray_shape ################
def check_ndarray_shape(v, vname, vshape, vaxis = None):
    generic_check_ndarray_shape(v, vname, vshape, vaxis)

def check_ndarray_axis_len(v, vname, vaxis, vlen):
    generic_check_ndarray_shape(v, vname, vlen, vaxis)

############# generic_check_ndarray_is_square ##############
def check_ndarray_is_square(v: Any, vname: str) -> None:
    generic_check_ndarray_is_square(v, vname)

############# check_ndarray_len_is_multiple_of #############
def check_ndarray_len_is_multiple_of(v, vname, m):
    if (len(v) % m) != 0:
        raise ValueError("len({0}) is not a multiple of {1}".format(vname, m))

def check_ndarray_len_is_multiple_of_2(v: Any, vname: str) -> None:
    check_ndarray_len_is_multiple_of(v, vname, 2)

def check_ndarray_len_is_multiple_of_3(v: Any, vname: str) -> None:
    check_ndarray_len_is_multiple_of(v, vname, 3)

def check_ndarray_len_is_multiple_of_4(v: Any, vname: str) -> None:
    check_ndarray_len_is_multiple_of(v, vname, 4)
