import numpy
from typing import Any
from typing import Union

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
def check_is_ndarray(v: Any, vname: str):
    generic_check_isinstance(v, vname, numpy.ndarray)

def check_is_Generator(v: Any, vname: str):
    generic_check_isinstance(v, vname, numpy.random.Generator)

def check_is_RandomState(v: Any, vname: str):
    generic_check_isinstance(v, vname, numpy.random.RandomState)

def check_is_Generator_or_RandomState(v: Any, vname: str):
    generic_check_isinstance(v, vname, (numpy.random.Generator,numpy.random.RandomState))

################################################################################
#################### conditional isinstance check functions ####################
################################################################################
def cond_check_is_ndarray(v, vname, cond = generic_default_cond):
    generic_cond_check_isinstance(v, vname, numpy.ndarray, cond)

def cond_check_is_Generator(v, vname, cond = generic_default_cond):
    generic_cond_check_isinstance(v, vname, numpy.random.Generator, cond)

def cond_check_is_RandomState(v, vname, cond = generic_default_cond):
    generic_cond_check_isinstance(v, vname, numpy.random.RandomState, cond)

def cond_check_is_Generator_or_RandomState(v, vname, cond = generic_default_cond):
    generic_cond_check_isinstance(v, vname, (numpy.random.Generator,numpy.random.RandomState), cond)

################################################################################
############################ dtype check functions #############################
################################################################################
def check_dtype(v: Any, vname: str, vdtype: Union[tuple,numpy.dtype]):
    istuple = isinstance(vdtype, tuple)
    logic = any(numpy.issubdtype(v, e) for e in vdtype) if istuple else numpy.issubdtype(v, vdtype)
    if not logic:
        raise TypeError("variable '{0}' must be of dtype '{1}'".format(vname, vdtype))

def check_dtype_is_bool(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("bool_")):
        raise TypeError("variable '{0}' must be of dtype 'bool'".format(vname))

def check_dtype_is_float16(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("float16")):
        raise TypeError("variable '{0}' must be of dtype 'float16'".format(vname))

def check_dtype_is_float32(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("float32")):
        raise TypeError("variable '{0}' must be of dtype 'float32'".format(vname))

def check_dtype_is_float64(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("float64")):
        raise TypeError("variable '{0}' must be of dtype 'float64'".format(vname))

def check_dtype_is_float128(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("float128")):
        raise TypeError("variable '{0}' must be of dtype 'float128'".format(vname))

def check_dtype_is_floating(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.floating):
        raise TypeError("variable '{0}' must be a floating dtype".format(vname))

def check_dtype_is_int8(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("int8")):
        raise TypeError("variable '{0}' must be of dtype 'int8'".format(vname))

def check_dtype_is_int16(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("int16")):
        raise TypeError("variable '{0}' must be of dtype 'int16'".format(vname))

def check_dtype_is_int32(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("int32")):
        raise TypeError("variable '{0}' must be of dtype 'int32'".format(vname))

def check_dtype_is_int64(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("int64")):
        raise TypeError("variable '{0}' must be of dtype 'int64'".format(vname))

def check_dtype_is_integer(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.integer):
        raise TypeError("variable '{0}' must be an integer dtype".format(vname))

def check_dtype_is_number(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.number):
        raise TypeError("variable '{0}' must be a numeric dtype".format(vname))

def check_dtype_is_object(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("object_")):
        raise TypeError("variable '{0}' must be of dtype 'object_'".format(vname))

def check_dtype_is_string(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("string_")):
        raise TypeError("variable '{0}' must be of dtype 'string_'".format(vname))

def check_dtype_is_unicode(v: Any, vname: str):
    if not numpy.issubdtype(v, numpy.dtype("unicode_")):
        raise TypeError("variable '{0}' must be of dtype 'unicode_'".format(vname))

################################################################################
######################## compound dtype check functions ########################
################################################################################
def check_dtype_is_bool_or_number(v: Any, vname: str):
    generic_check_dtype_issubdtype(v, vname, (numpy.bool_, numpy.number))

def check_dtype_is_integer_or_floating(v: Any, vname: str):
    generic_check_dtype_issubdtype(v, vname, (numpy.integer, numpy.floating))

def check_dtype_is_object_or_string(v: Any, vname: str):
    generic_check_dtype_issubdtype(v, vname, (numpy.object_, numpy.string_))

################################################################################
######################## ndarray dtype check functions #########################
################################################################################
def check_ndarray_dtype(v, vname, vdtype):
    generic_check_ndarray_dtype_issubdtype(v, vname, vdtype)

def check_ndarray_dtype_is_bool(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.bool_)

def check_ndarray_dtype_is_float16(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.float16)

def check_ndarray_dtype_is_float32(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.float32)

def check_ndarray_dtype_is_float64(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.float64)

def check_ndarray_dtype_is_float128(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.float128)

def check_ndarray_dtype_is_floating(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.floating)

def check_ndarray_dtype_is_int8(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.int8)

def check_ndarray_dtype_is_int16(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.int16)

def check_ndarray_dtype_is_int32(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.int32)

def check_ndarray_dtype_is_int64(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.int64)

def check_ndarray_dtype_is_integer(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.integer)

def check_ndarray_dtype_is_number(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.number)

def check_ndarray_dtype_is_object(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.object_)

def check_ndarray_dtype_is_string(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.string_)

def check_ndarray_dtype_is_unicode(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, numpy.unicode_)

################################################################################
#################### compound ndarray dtype check functions ####################
################################################################################
def check_ndarray_dtype_is_bool_or_number(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, (numpy.bool_, numpy.number))

def check_ndarray_dtype_is_integer_or_floating(v: Any, vname: str):
    generic_check_ndarray_dtype_issubdtype(v, vname, (numpy.integer, numpy.floating))

def check_ndarray_dtype_is_object_or_string(v: Any, vname: str):
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
