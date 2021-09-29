import numpy

from . import generic_check_ndarray_eq
from . import generic_check_ndarray_sum
from . import generic_check_ndarray_ndim
from . import generic_check_ndarray_size
from . import generic_check_ndarray_shape
from . import generic_check_ndarray_is_square

from . import generic_default_cond
from . import generic_cond_check_ndarray_eq
from . import generic_cond_check_ndarray_sum
from . import generic_cond_check_ndarray_ndim
from . import generic_cond_check_ndarray_size
from . import generic_cond_check_ndarray_shape
from . import generic_cond_check_ndarray_is_square
from . import generic_check_ndarray_ndim_gteq

################################################################################
############################### check functions ################################
################################################################################

################# generic_check_ndarray_eq #################
def check_ndarray_eq(v, vname, w, wname):
    generic_check_ndarray_eq(v, vname, w, wname)

def check_ndarray_is_binary(v, vname):
    if not numpy.all((v == 0) | (v == 1)):
        raise ValueError("variable '{0}' must have values equal to 0 or 1".format(vname))

def check_ndarray_is_positive(v, vname):
    if numpy.any(v < 0):
        raise ValueError("variable '{0}' must have all positive values".format(vname))

################ generic_check_ndarray_ndim ################
def check_ndarray_ndim(v, vname, vndim):
    generic_check_ndarray_ndim(v, vname, vndim)

def check_ndarray_ndim_gteq(v, vname, vndim):
    generic_check_ndarray_ndim_gteq(v, vname, vndim)

def check_ndarray_is_1d(v, vname):
    generic_check_ndarray_ndim(v, vname, 1)

def check_ndarray_at_least_1d(v, vname):
    generic_check_ndarray_ndim_gteq(v, vname, 1)

def check_ndarray_is_2d(v, vname):
    generic_check_ndarray_ndim(v, vname, 2)

def check_ndarray_at_least_2d(v, vname):
    generic_check_ndarray_ndim_gteq(v, vname, 2)

def check_ndarray_is_3d(v, vname):
    generic_check_ndarray_ndim(v, vname, 3)

def check_ndarray_at_least_3d(v, vname):
    generic_check_ndarray_ndim_gteq(v, vname, 3)

################ generic_check_ndarray_size ################
def check_ndarray_size(v, vname, vsize):
    generic_check_ndarray_size(v, vname, vsize)

################ generic_check_ndarray_sum #################
def check_ndarray_sum(v, vname, vsum, vaxis):
    generic_check_ndarray_sum(v, vname, vsum, vaxis)

############### generic_check_ndarray_shape ################
def check_ndarray_shape(v, vname, vshape, vaxis = None):
    generic_check_ndarray_shape(v, vname, vshape, vaxis)

def check_ndarray_axis_len(v, vname, vaxis, vlen):
    generic_check_ndarray_shape(v, vname, vlen, vaxis)

############# generic_check_ndarray_is_square ##############
def check_ndarray_is_square(v, vname):
    generic_check_ndarray_is_square(v, vname)

############# check_ndarray_len_is_multiple_of #############
def check_ndarray_len_is_multiple_of(v, vname, m):
    if (len(v) % m) != 0:
        raise ValueError("len({0}) is not a multiple of {1}".format(vname, m))

def check_ndarray_len_is_multiple_of_2(v, vname):
    check_ndarray_len_is_multiple_of(v, vname, 2)

def check_ndarray_len_is_multiple_of_3(v, vname):
    check_ndarray_len_is_multiple_of(v, vname, 3)

def check_ndarray_len_is_multiple_of_4(v, vname):
    check_ndarray_len_is_multiple_of(v, vname, 4)

################################################################################
######################### conditional check functions ##########################
################################################################################

################# generic_check_ndarray_eq #################
def cond_check_ndarray_eq(v, vname, w, wname, cond = generic_default_cond):
    generic_cond_check_ndarray_eq(v, vname, w, wname, cond)

def cond_check_ndarray_is_binary(v, vname, cond = generic_default_cond):
    if cond(v):
        check_ndarray_is_binary(v, vname)

################ generic_check_ndarray_ndim ################
def cond_check_ndarray_ndim(v, vname, vndim, cond = generic_default_cond):
    generic_cond_check_ndarray_ndim(v, vname, vndim, cond)

def cond_check_ndarray_is_1d(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_ndim(v, vname, 1, cond)

def cond_check_ndarray_is_2d(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_ndim(v, vname, 2, cond)

def cond_check_ndarray_is_3d(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_ndim(v, vname, 3, cond)

################ generic_check_ndarray_size ################
def cond_check_ndarray_size(v, vname, vsize, cond = generic_default_cond):
    generic_cond_check_ndarray_size(v, vname, vsize, cond)

################ generic_check_ndarray_sum #################
def cond_check_ndarray_sum(v, vname, vsum, vaxis, cond = generic_default_cond):
    generic_cond_check_ndarray_sum(v, vname, vsum, vaxis, cond)

############# generic_cond_check_ndarray_shape #############
def cond_check_ndarray_shape(v, vname, vshape, vaxis = None, cond = generic_default_cond):
    generic_cond_check_ndarray_shape(v, vname, vshape, vaxis, cond)

def cond_check_ndarray_axis_len(v, vname, vaxis, vlen, cond = generic_default_cond):
    generic_cond_check_ndarray_shape(v, vname, vlen, vaxis, cond)

########### generic_cond_check_ndarray_is_square ###########
def check_ndarray_is_square(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_is_square(v, vname, cond)
