"""
Module containing subroutines to check Python types.
"""

__all__ = [
    "check_inherits",
    "check_isinstance",
    "check_is_bool",
    "check_is_class",
    "check_is_dict",
    "check_is_float",
    "check_is_int",
    "check_is_int_or_inf",
    "check_is_int_or_None",
    "check_is_list",
    "check_is_range",
    "check_is_set",
    "check_is_str",
    "check_is_tuple",
    "check_is_type",
    "check_is_Complex",
    "check_is_Integral",
    "check_is_Number",
    "check_is_Real",
    "check_is_Callable",
    "check_is_Container",
    "check_is_Sequence",
    "check_is_array_like",
    "check_is_str_or_iterable",
    "check_is_list_or_tuple",
]

import math
from numbers import Complex
from numbers import Integral
from numbers import Number
from numbers import Real
from typing import Callable
from typing import Container
from typing import Sequence
from typing import Union
import inspect

################################################################################
###################### basic inheritance check functions #######################
################################################################################
def check_inherits(v: object, vname: str, vtype: type) -> None:
    """
    Generic check of inheritance using method resolution order metadata.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the object variable.
    vtype : type
        Object type from which 'obj' must inherit.
    """
    if isinstance(vtype, type):
        if vtype not in v.__mro__:
            raise TypeError("variable '{0}' must inherit from '{1}'".format(vname, vtype))
    elif isinstance(vtype, tuple):
        if all(e not in v.__mro__ for e in vtype):
            raise TypeError("variable '{0}' must inherit from all of '{1}'".format(vname, vtype))
    else:
        raise TypeError("'objtype' must be of type 'type' or 'tuple'")

################################################################################
################## basic check functions for basic data types ##################
################################################################################
def check_isinstance(v: object, vname: str, vtype: Union[type,tuple]) -> None:
    """
    Generic check type function.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    vtype : type, tuple
        type : Type associated with the object variable.
        tuple : Tuple of acceptable types (logical or) for the object variable.
    """
    if not isinstance(v, vtype):
        tname = None
        if isinstance(vtype, tuple):
            tname = ""
            l = len(vtype)
            for e, i in enumerate(vtype):
                tname += "'" + e.__name__ + "'"
                if i < l - 2:
                    tname += ", "
                elif i < l - 1:
                    tname += ", or "
        else:
            tname = vtype.__name__
        raise TypeError("variable '{0}' must of type {1}".format(vname, tname))

def check_is_bool(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``bool``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, bool):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,bool.__name__,type(v).__name__))

def check_is_bytes(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``bytes``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, bytes):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,bytes.__name__,type(v).__name__))

def check_is_class(v: object, vname: str) -> None:
    """
    Check if a Python object is a ``class``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not inspect.isclass(v):
        raise TypeError("variable '{0}' must be a class name".format(vname))

def check_is_dict(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``dict``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, dict):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,dict.__name__,type(v).__name__))

def check_is_float(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``float``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, float):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,float.__name__,type(v).__name__))

def check_is_frozenset(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``frozenset``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, frozenset):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,frozenset.__name__,type(v).__name__))

def check_is_int(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``int``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, int):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,int.__name__,type(v).__name__))

def check_is_list(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``list``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, list):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,list.__name__,type(v).__name__))

def check_is_range(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``range``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, range):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,range.__name__,type(v).__name__))

def check_is_set(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``set``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, set):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,set.__name__,type(v).__name__))

def check_is_str(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``str``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, str):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,str.__name__,type(v).__name__))

def check_is_tuple(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``tuple``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, tuple):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,tuple.__name__,type(v).__name__))

def check_is_type(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``type``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, type):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,type.__name__,type(v).__name__))

################################################################################
#################### basic check functions for number types ####################
################################################################################
def check_is_Complex(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``numbers.Complex``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, Complex):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,Complex.__name__,type(v).__name__))

def check_is_Integral(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``numbers.Integral``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, Integral):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,Integral.__name__,type(v).__name__))

def check_is_Number(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``numbers.Number``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, Number):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,Number.__name__,type(v).__name__))

def check_is_Real(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``numbers.Real``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, Real):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,Real.__name__,type(v).__name__))

################################################################################
#################### basic check functions for typing types ####################
################################################################################
def check_is_Callable(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``typing.Callable``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, Callable):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,Callable.__name__,type(v).__name__))

def check_is_Container(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``typing.Container``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, Container):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,Container.__name__,type(v).__name__))

def check_is_Sequence(v: object, vname: str) -> None:
    """
    Check if a Python object is of type ``typing.Sequence``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, Sequence):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,Sequence.__name__,type(v).__name__))

################################################################################
################ compound check functions for basic data types #################
################################################################################
def check_is_array_like(v: object, vname: str) -> None:
    """
    Check if a Python object is array-like (has ``__len__``, ``__iter__``, and 
    ``__getitem__`` atributes). Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    alattr = ("__len__","__iter__","__getitem__")
    for a in alattr:
        if not hasattr(v, a):
            raise TypeError("'{0}' must have attribute '{1}' to be array_like".format(vname,a))

def check_is_str_or_iterable(v: object, vname: str) -> None:
    """
    Check if a Python object is a ``str`` or is iterable. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not (isinstance(v, str) or hasattr(v, "__iter__")):
        raise TypeError("'{0}' must be of type str or have attribute __iter__.".format(vname))

def check_is_list_or_tuple(v: object, vname: str) -> None:
    """
    Check if a Python object is a ``list`` or ``tuple``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, (list,tuple)):
        raise TypeError("variable '{0}' must be of type '{1}' or '{2}' but received type '{3}'".format(vname,list.__name__,tuple.__name__,type(v).__name__))

def check_is_str_or_Integral(v: object, vname: str) -> None:
    """
    Check if a Python object is a ``str`` or ``Integral``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, (str,Integral)):
        raise TypeError(
            "variable '{0}' must be of type '{1}' or '{2}' but received type '{3}'".format(
                vname,
                str.__name__,
                Integral.__name__,
                type(v).__name__
            )
        )

def check_is_str_or_Sequence(v: object, vname: str) -> None:
    """
    Check if a Python object is a ``str`` or ``Sequence``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, (str,Sequence)):
        raise TypeError(
            "variable '{0}' must be of type '{1}' or '{2}' but received type '{3}'".format(
                vname,
                str.__name__,
                Sequence.__name__,
                type(v).__name__
            )
        )

def check_is_int_or_inf(v: object, vname: str) -> None:
    """
    Check if a Python object is an ``int`` or if it is infinite. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if (v != math.inf) and (not isinstance(v, int)):
        raise TypeError("variable '{0}' must be of type 'int' or be infinite".format(vname))

def check_is_int_or_None(v: object, vname: str) -> None:
    """
    Check if a Python object is an ``int`` or if it is ``None``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if (not isinstance(v, int)) and (v is not None):
        raise TypeError("variable '{0}' must be of type 'int' or None".format(vname))

def check_is_Integral_or_inf(v: object, vname: str) -> None:
    """
    Check if a Python object is a ``numbers.Integral`` or infinite. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if (v != math.inf) and (not isinstance(v, Integral)):
        raise TypeError("variable '{0}' must be of type '{1}' or the value 'math.inf' but received type '{2}'".format(vname,Integral.__name__,type(v).__name__))

def check_is_Integral_or_None(v: object, vname: str) -> None:
    """
    Check if a Python object is a ``numbers.Integral`` or ``None``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if (v is not None) and (not isinstance(v, Integral)):
        raise TypeError("variable '{0}' must be of type '{1}' or the value 'None' but received type '{2}'".format(vname,Integral.__name__,type(v).__name__))

### Iterating type checks ###

def check_Sequence_all_type(v: Sequence, vname: str, vtype: Union[type,tuple]) -> None:
    if not all(isinstance(e, vtype) for e in v):
        raise TypeError(
            "Sequence '{0}' must have values all of type '{1}'".format(
                vname,
                vtype.__name__
            )
        )

