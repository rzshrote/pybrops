import numpy

from . import generic_check_dtype_issubdtype
from . import generic_check_isinstance
from . import generic_check_ndarray_dtype_issubdtype

from . import generic_default_cond
from . import generic_cond_check_dtype_issubdtype
from . import generic_cond_check_isinstance
from . import generic_cond_check_ndarray_dtype_issubdtype

################################################################################
########################## isinstance check functions ##########################
################################################################################
def check_is_ndarray(v, vname):
    generic_check_isinstance(v, vname, numpy.ndarray)

def check_is_Generator(v, vname):
    generic_check_isinstance(v, vname, numpy.random.Generator)

################################################################################
#################### conditional isinstance check functions ####################
################################################################################
def cond_check_is_ndarray(v, vname, cond = generic_default_cond):
    generic_cond_check_isinstance(v, vname, numpy.ndarray, cond)

def check_is_Generator(v, vname, cond = generic_default_cond):
    generic_cond_check_isinstance(v, vname, numpy.random.Generator, cond)

################################################################################
############################ dtype check functions #############################
################################################################################
def check_dtype(v, vname, vdtype):
    generic_check_dtype_issubdtype(v, vname, vdtype)

def check_dtype_is_bool(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.bool_)

def check_dtype_is_float16(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.float16)

def check_dtype_is_float32(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.float32)

def check_dtype_is_float64(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.float64)

def check_dtype_is_float128(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.float128)

def check_dtype_is_floating(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.floating)

def check_dtype_is_int8(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.int8)

def check_dtype_is_int16(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.int16)

def check_dtype_is_int32(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.int32)

def check_dtype_is_int64(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.int64)

def check_dtype_is_integer(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.integer)

def check_dtype_is_number(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.number)

def check_dtype_is_object(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.object_)

def check_dtype_is_string(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.string_)

def check_dtype_is_unicode(v, vname):
    generic_check_dtype_issubdtype(v, vname, numpy.unicode_)

################################################################################
######################## compound dtype check functions ########################
################################################################################
def check_dtype_is_bool_or_number(v, vname):
    generic_check_dtype_issubdtype(v, vname, (numpy.bool_, numpy.number))

def check_dtype_is_integer_or_floating(v, vname):
    generic_check_dtype_issubdtype(v, vname, (numpy.integer, numpy.floating))

def check_dtype_is_object_or_string(v, vname):
    generic_check_dtype_issubdtype(v, vname, (numpy.object_, numpy.string_))

################################################################################
###################### conditional dtype check functions #######################
################################################################################
def cond_check_dtype(v, vname, vdtype, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, vdtype, cond)

def cond_check_dtype_is_bool(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.bool_, cond)

def cond_check_dtype_is_floating(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.floating, cond)

def cond_check_dtype_is_float16(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.float16, cond)

def cond_check_dtype_is_float32(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.float32, cond)

def cond_check_dtype_is_float64(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.float64, cond)

def cond_check_dtype_is_float128(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.float128, cond)

def cond_check_dtype_is_floating(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.floating, cond)

def cond_check_dtype_is_int8(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.int8, cond)

def cond_check_dtype_is_int16(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.int16, cond)

def cond_check_dtype_is_int32(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.int32, cond)

def cond_check_dtype_is_int64(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.int64, cond)

def cond_check_dtype_is_integer(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.integer, cond)

def cond_check_dtype_is_number(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.number, cond)

def cond_check_dtype_is_object(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.object_, cond)

def cond_check_dtype_is_string(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.string_, cond)

def cond_check_dtype_is_unicode(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, numpy.unicode_, cond)

################################################################################
################## conditional compound dtype check functions ##################
################################################################################
def cond_check_dtype_is_bool_or_number(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, (numpy.bool_, numpy.number), cond)

def cond_check_dtype_is_integer_or_floating(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, (numpy.integer, numpy.floating), cond)

def cond_check_dtype_is_object_or_string(v, vname, cond = generic_default_cond):
    generic_cond_check_dtype_issubdtype(v, vname, (numpy.object_, numpy.string_), cond)

################################################################################
######################## ndarray dtype check functions #########################
################################################################################
def check_ndarray_dtype(v, vname, vdtype):
    generic_check_ndarray_dtype_issubdtype(v, vname, vdtype)

def check_ndarray_dtype_is_bool(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.bool_)

def check_ndarray_dtype_is_float16(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.float16)

def check_ndarray_dtype_is_float32(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.float32)

def check_ndarray_dtype_is_float64(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.float64)

def check_ndarray_dtype_is_float128(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.float128)

def check_ndarray_dtype_is_floating(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.floating)

def check_ndarray_dtype_is_int8(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.int8)

def check_ndarray_dtype_is_int16(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.int16)

def check_ndarray_dtype_is_int32(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.int32)

def check_ndarray_dtype_is_int64(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.int64)

def check_ndarray_dtype_is_integer(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.integer)

def check_ndarray_dtype_is_number(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.number)

def check_ndarray_dtype_is_object(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.object_)

def check_ndarray_dtype_is_string(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.string_)

def check_ndarray_dtype_is_unicode(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.unicode_)

################################################################################
#################### compound ndarray dtype check functions ####################
################################################################################
def check_ndarray_dtype_is_bool_or_number(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, (numpy.bool_, numpy.number))

def check_ndarray_dtype_is_integer_or_floating(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, (numpy.integer, numpy.floating))

def check_ndarray_dtype_is_object_or_string(v, vname):
    generic_check_ndarray_dtype_issubdtype(v, vname, (numpy.object_, numpy.string_))

################################################################################
################## conditional ndarray dtype check functions ###################
################################################################################
def cond_check_ndarray_dtype(v, vname, vdtype, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, vdtype, cond)

def cond_check_ndarray_dtype_is_bool(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.bool_, cond)

def cond_check_ndarray_dtype_is_floating(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.floating, cond)

def cond_check_ndarray_dtype_is_float16(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.float16, cond)

def cond_check_ndarray_dtype_is_float32(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.float32, cond)

def cond_check_ndarray_dtype_is_float64(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.float64, cond)

def cond_check_ndarray_dtype_is_float128(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.float128, cond)

def cond_check_ndarray_dtype_is_floating(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.floating, cond)

def cond_check_ndarray_dtype_is_int8(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.int8, cond)

def cond_check_ndarray_dtype_is_int16(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.int16, cond)

def cond_check_ndarray_dtype_is_int32(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.int32, cond)

def cond_check_ndarray_dtype_is_int64(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.int64, cond)

def cond_check_ndarray_dtype_is_integer(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.integer, cond)

def cond_check_ndarray_dtype_is_number(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.number, cond)

def cond_check_ndarray_dtype_is_object(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.object_, cond)

def cond_check_ndarray_dtype_is_string(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.string_, cond)

def cond_check_ndarray_dtype_is_unicode(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, numpy.unicode_, cond)

################################################################################
############## conditional compound ndarray dtype check functions ##############
################################################################################
def cond_check_ndarray_dtype_is_bool_or_number(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, (numpy.bool_, numpy.number), cond)

def cond_check_ndarray_dtype_is_integer_or_floating(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, (numpy.integer, numpy.floating), cond)

def cond_check_ndarray_dtype_is_object_or_string(v, vname, cond = generic_default_cond):
    generic_cond_check_ndarray_dtype_issubdtype(v, vname, (numpy.object_, numpy.string_), cond)
