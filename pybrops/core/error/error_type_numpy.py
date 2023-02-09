import numpy
from typing import Any
from typing import Union

################################################################################
########################## isinstance check functions ##########################
################################################################################
def check_is_ndarray(v: Any, vname: str) -> None:
    if not isinstance(v, numpy.ndarray):
        raise TypeError("variable '{0}' must be of type 'numpy.ndarray'".format(vname))

def check_is_Generator(v: Any, vname: str) -> None:
    if not isinstance(v, numpy.random.Generator):
        raise TypeError("variable '{0}' must be of type 'numpy.random.Generator'".format(vname))

def check_is_RandomState(v: Any, vname: str) -> None:
    if not isinstance(v, numpy.random.RandomState):
        raise TypeError("variable '{0}' must be of type 'numpy.random.RandomState'".format(vname))

def check_is_Generator_or_RandomState(v: Any, vname: str) -> None:
    if not (isinstance(v, numpy.random.Generator) or isinstance(v, numpy.random.RandomState)):
        raise TypeError("variable '{0}' must be of type 'numpy.random.Generator' or 'numpy.random.RandomState'".format(vname))

################################################################################
############################ dtype check functions #############################
################################################################################
def check_dtype(v: Union[numpy.dtype,str], vname: str, vdtype: Union[tuple,numpy.dtype]):
    istuple = isinstance(vdtype, tuple)
    logic = any(numpy.issubdtype(v, e) for e in vdtype) if istuple else numpy.issubdtype(v, vdtype)
    if not logic:
        raise TypeError("variable '{0}' must be of dtype '{1}'".format(vname, vdtype))

def check_dtype_is_bool(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("bool_")):
        raise TypeError("variable '{0}' must be of dtype 'bool'".format(vname))

def check_dtype_is_float16(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("float16")):
        raise TypeError("variable '{0}' must be of dtype 'float16'".format(vname))

def check_dtype_is_float32(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("float32")):
        raise TypeError("variable '{0}' must be of dtype 'float32'".format(vname))

def check_dtype_is_float64(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("float64")):
        raise TypeError("variable '{0}' must be of dtype 'float64'".format(vname))

def check_dtype_is_float128(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("float128")):
        raise TypeError("variable '{0}' must be of dtype 'float128'".format(vname))

def check_dtype_is_floating(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.floating):
        raise TypeError("variable '{0}' must be a floating dtype".format(vname))

def check_dtype_is_int8(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("int8")):
        raise TypeError("variable '{0}' must be of dtype 'int8'".format(vname))

def check_dtype_is_int16(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("int16")):
        raise TypeError("variable '{0}' must be of dtype 'int16'".format(vname))

def check_dtype_is_int32(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("int32")):
        raise TypeError("variable '{0}' must be of dtype 'int32'".format(vname))

def check_dtype_is_int64(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("int64")):
        raise TypeError("variable '{0}' must be of dtype 'int64'".format(vname))

def check_dtype_is_integer(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.integer):
        raise TypeError("variable '{0}' must be an integer dtype".format(vname))

def check_dtype_is_number(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.number):
        raise TypeError("variable '{0}' must be a numeric dtype".format(vname))

def check_dtype_is_object(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("object_")):
        raise TypeError("variable '{0}' must be of dtype 'object_'".format(vname))

def check_dtype_is_string(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("string_")):
        raise TypeError("variable '{0}' must be of dtype 'string_'".format(vname))

def check_dtype_is_unicode(v: Union[numpy.dtype,str], vname: str):
    if not numpy.issubdtype(v, numpy.dtype("unicode_")):
        raise TypeError("variable '{0}' must be of dtype 'unicode_'".format(vname))

################################################################################
######################## compound dtype check functions ########################
################################################################################
def check_dtype_is_bool_or_number(v: Union[numpy.dtype,str], vname: str):
    if not (numpy.issubdtype(v, numpy.dtype("bool_")) or numpy.issubdtype(v, numpy.number)):
        raise TypeError("variable '{0}' must be of dtype 'bool' or be numeric".format(vname))

def check_dtype_is_integer_or_floating(v: Union[numpy.dtype,str], vname: str):
    if not (numpy.issubdtype(v, numpy.integer) or numpy.issubdtype(v, numpy.floating)):
        raise TypeError("variable '{0}' must be an integer or floating dtype".format(vname))

def check_dtype_is_object_or_string(v: Union[numpy.dtype,str], vname: str):
    if not (numpy.issubdtype(v, numpy.dtype("object_")) or numpy.issubdtype(v, numpy.dtype("string_"))):
        raise TypeError("variable '{0}' must be of dtype 'object_' or 'string_'".format(vname))

################################################################################
######################## ndarray dtype check functions #########################
################################################################################
def check_ndarray_dtype(v: numpy.ndarray, vname: str, vdtype: Union[tuple,numpy.dtype,str]):
    istuple = isinstance(vdtype, tuple)
    logic = any(numpy.issubdtype(v.dtype, e) for e in vdtype) if istuple else numpy.issubdtype(v.dtype, vdtype)
    if not logic:
        raise TypeError("numpy.ndarray '{0}' must be of dtype '{1}'".format(vname, vdtype))

def check_ndarray_dtype_is_bool(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("bool_")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'bool'".format(vname))

def check_ndarray_dtype_is_float16(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("float16")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'float16'".format(vname))

def check_ndarray_dtype_is_float32(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("float32")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'float32'".format(vname))

def check_ndarray_dtype_is_float64(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("float64")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'float64'".format(vname))

def check_ndarray_dtype_is_float128(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("float128")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'float128'".format(vname))

def check_ndarray_dtype_is_floating(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.floating):
        raise TypeError("numpy.ndarray '{0}' must be a floating dtype".format(vname))

def check_ndarray_dtype_is_int8(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("int8")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'int8'".format(vname))

def check_ndarray_dtype_is_int16(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("int16")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'int16'".format(vname))

def check_ndarray_dtype_is_int32(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("int32")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'int32'".format(vname))

def check_ndarray_dtype_is_int64(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("int64")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'int64'".format(vname))

def check_ndarray_dtype_is_integer(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.integer):
        raise TypeError("numpy.ndarray '{0}' must be an integer dtype".format(vname))

def check_ndarray_dtype_is_number(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.number):
        raise TypeError("numpy.ndarray '{0}' must be a numeric dtype".format(vname))

def check_ndarray_dtype_is_object(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("object_")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'object_'".format(vname))

def check_ndarray_dtype_is_string(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("string_")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'string_'".format(vname))

def check_ndarray_dtype_is_unicode(v: numpy.ndarray, vname: str):
    if not numpy.issubdtype(v.dtype, numpy.dtype("unicode_")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'unicode_'".format(vname))

################################################################################
#################### compound ndarray dtype check functions ####################
################################################################################
def check_ndarray_dtype_is_bool_or_number(v: numpy.ndarray, vname: str):
    if not (numpy.issubdtype(v.dtype, numpy.dtype("bool_")) or numpy.issubdtype(v.dtype, numpy.number)):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'bool' or be numeric".format(vname))

def check_ndarray_dtype_is_integer_or_floating(v: numpy.ndarray, vname: str):
    if not (numpy.issubdtype(v.dtype, numpy.integer) or numpy.issubdtype(v.dtype, numpy.floating)):
        raise TypeError("numpy.ndarray '{0}' must be an integer or floating dtype".format(vname))

def check_ndarray_dtype_is_object_or_string(v: numpy.ndarray, vname: str):
    if not (numpy.issubdtype(v.dtype, numpy.dtype("object_")) or numpy.issubdtype(v.dtype, numpy.dtype("string_"))):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'object_' or 'string_'".format(vname))
