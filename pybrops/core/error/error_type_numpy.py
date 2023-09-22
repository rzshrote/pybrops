"""
Module containing subroutines to check ``numpy`` types.
"""

__all__ = [
    "check_is_ndarray",
    "check_is_Generator",
    "check_is_RandomState",
    "check_is_Generator_or_RandomState",
    "check_is_Integral_or_ndarray",
    "check_is_Number_or_ndarray",
    "check_is_Real_or_ndarray",
    "check_ndarray_dtype",
    "check_ndarray_dtype_is_bool",
    "check_ndarray_dtype_is_float16",
    "check_ndarray_dtype_is_float32",
    "check_ndarray_dtype_is_float64",
    "check_ndarray_dtype_is_floating",
    "check_ndarray_dtype_is_int8",
    "check_ndarray_dtype_is_int16",
    "check_ndarray_dtype_is_int32",
    "check_ndarray_dtype_is_int64",
    "check_ndarray_dtype_is_integer",
    "check_ndarray_dtype_is_number",
    "check_ndarray_dtype_is_object",
    "check_ndarray_dtype_is_real",
    "check_ndarray_dtype_is_string",
    "check_ndarray_dtype_is_unicode",
    "check_ndarray_dtype_is_bool_or_integer",
    "check_ndarray_dtype_is_bool_or_number",
    "check_ndarray_dtype_is_integer_or_floating",
    "check_ndarray_dtype_is_object_or_string",
]

from numbers import Integral, Number, Real
import numpy
from numpy.random import Generator, RandomState
from typing import Union

################################################################################
########################## isinstance check functions ##########################
################################################################################
def check_is_str_or_ndarray(v: object, vname: str) -> None:
    """
    Check whether a Python object is a ``str`` or ``numpy.ndarray``.

    Parameters
    ----------
    v : object
        Python object to test.
    vname : str
        Name assigned to the Python object for an error message.
    """
    if not isinstance(v, (str,numpy.ndarray)):
        raise TypeError(
            "variable '{0}' must be of type '{1}' or '{2}' but received type '{3}'".format(
                vname,
                str.__name__,
                numpy.ndarray.__name__,
                type(v).__name__
            )
        )

def check_is_ndarray(v: object, vname: str) -> None:
    """
    Check whether a Python object is a ``numpy.ndarray``.

    Parameters
    ----------
    v : object
        Python object to test.
    vname : str
        Name assigned to the Python object for an error message.
    """
    if not isinstance(v, numpy.ndarray):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,numpy.ndarray.__name__,type(v).__name__))

def check_is_Generator(v: object, vname: str) -> None:
    """
    Check whether a Python object is a ``numpy.random.Generator``.

    Parameters
    ----------
    v : object
        Python object to test.
    vname : str
        Name assigned to the Python object for an error message.
    """
    if not isinstance(v, Generator):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,Generator.__name__,type(v).__name__))

def check_is_RandomState(v: object, vname: str) -> None:
    """
    Check whether a Python object is a ``numpy.random.RandomState``.

    Parameters
    ----------
    v : object
        Python object to test.
    vname : str
        Name assigned to the Python object for an error message.
    """
    if not isinstance(v, RandomState):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,RandomState.__name__,type(v).__name__))

def check_is_Generator_or_RandomState(v: object, vname: str) -> None:
    """
    Check whether a Python object is a ``numpy.random.Generator`` or ``numpy.random.RandomState``.

    Parameters
    ----------
    v : object
        Python object to test.
    vname : str
        Name assigned to the Python object for an error message.
    """
    if not (isinstance(v, Generator) or isinstance(v, RandomState)):
        raise TypeError("variable '{0}' must be of type '{1}' or '{2}' but received type '{3}'".format(vname,Generator.__name__,RandomState.__name__,type(v).__name__))

################################################################################
##################### compound isinstance check functions ######################
################################################################################
def check_is_Integral_or_ndarray(v: object, vname: str) -> None:
    """
    Check whether a Python object is a ``numpy.random.Generator`` or ``numpy.random.RandomState``.

    Parameters
    ----------
    v : object
        Python object to test.
    vname : str
        Name assigned to the Python object for an error message.
    """
    if not isinstance(v, (Integral,numpy.ndarray)):
        raise TypeError("variable '{0}' must be of type '{1}' or '{2}' but received type '{3}'".format(vname,Integral.__name__,numpy.ndarray.__name__,type(v).__name__))

def check_is_Number_or_ndarray(v: object, vname: str) -> None:
    """
    Check whether a Python object is a ``numpy.random.Generator`` or ``numpy.random.RandomState``.

    Parameters
    ----------
    v : object
        Python object to test.
    vname : str
        Name assigned to the Python object for an error message.
    """
    if not isinstance(v, (Number,numpy.ndarray)):
        raise TypeError("variable '{0}' must be of type '{1}' or '{2}' but received type '{3}'".format(vname,Number.__name__,numpy.ndarray.__name__,type(v).__name__))

def check_is_Real_or_ndarray(v: object, vname: str) -> None:
    """
    Check whether a Python object is a ``numpy.random.Generator`` or ``numpy.random.RandomState``.

    Parameters
    ----------
    v : object
        Python object to test.
    vname : str
        Name assigned to the Python object for an error message.
    """
    if not isinstance(v, (Real,numpy.ndarray)):
        raise TypeError("variable '{0}' must be of type '{1}' or '{2}' but received type '{3}'".format(vname,Real.__name__,numpy.ndarray.__name__,type(v).__name__))

################################################################################
######################## ndarray dtype check functions #########################
################################################################################
def check_ndarray_dtype(v: numpy.ndarray, vname: str, vdtype: Union[tuple,numpy.dtype,str]):
    istuple = isinstance(vdtype, tuple)
    logic = any(numpy.issubdtype(v.dtype, e) for e in vdtype) if istuple else numpy.issubdtype(v.dtype, vdtype)
    if not logic:
        raise TypeError("numpy.ndarray '{0}' must be of dtype '{1}'".format(vname, vdtype))

def check_ndarray_dtype_is_bool(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``bool`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.dtype("bool_")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'bool' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_float16(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``float16`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.dtype("float16")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'float16' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_float32(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``float32`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.dtype("float32")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'float32' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_float64(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``float64`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.dtype("float64")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'float64' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_floating(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``floating`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.floating):
        raise TypeError("numpy.ndarray '{0}' must be a floating dtype but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_int8(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``int8`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.dtype("int8")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'int8' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_int16(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``int16`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.dtype("int16")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'int16' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_int32(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``int32`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.dtype("int32")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'int32' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_int64(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``int64`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.dtype("int64")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'int64' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_integer(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``integer`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.integer):
        raise TypeError("numpy.ndarray '{0}' must be an integer dtype but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_number(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``number`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.number):
        raise TypeError("numpy.ndarray '{0}' must be a numeric dtype but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_object(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``object`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.dtype("object_")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'object_' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_real(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``floating`` or ``integer`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not (numpy.issubdtype(v.dtype, numpy.floating) or numpy.issubdtype(v.dtype, numpy.integer)):
        raise TypeError("numpy.ndarray '{0}' must be a real (floating or integer) dtype but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_string(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``string_`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.dtype("string_")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'string_' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_unicode(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``unicode_`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not numpy.issubdtype(v.dtype, numpy.dtype("unicode_")):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'unicode_' but received dtype '{1}'".format(vname,v.dtype))

################################################################################
#################### compound ndarray dtype check functions ####################
################################################################################
def check_ndarray_dtype_is_bool_or_integer(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``bool`` or ``integer`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not (numpy.issubdtype(v.dtype, numpy.dtype(bool)) or numpy.issubdtype(v.dtype, numpy.integer)):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'bool' or 'integer' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_bool_or_number(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``bool`` or ``number`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not (numpy.issubdtype(v.dtype, numpy.dtype(bool)) or numpy.issubdtype(v.dtype, numpy.number)):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'bool' or 'number' but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_integer_or_floating(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``integer`` or ``floating`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not (numpy.issubdtype(v.dtype, numpy.integer) or numpy.issubdtype(v.dtype, numpy.floating)):
        raise TypeError("numpy.ndarray '{0}' must be an integer or floating dtype but received dtype '{1}'".format(vname,v.dtype))

def check_ndarray_dtype_is_object_or_string(v: numpy.ndarray, vname: str) -> None:
    """
    Check if a ``numpy.ndarray`` has a ``object`` or ``string_`` dtype.

    Parameters
    ----------
    v : numpy.ndarray
        Array for which to check the ``dtype``.
    vname : str
        Name of the ``numpy.ndarray`` to use in the error message.
    """
    if not (numpy.issubdtype(v.dtype, numpy.dtype(object)) or numpy.issubdtype(v.dtype, numpy.dtype("string_"))):
        raise TypeError("numpy.ndarray '{0}' must be of dtype 'object' or 'string_' but received dtype '{1}'".format(vname,v.dtype))
