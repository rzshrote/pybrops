import math
from numbers import Complex, Integral, Number, Real
from typing import Any, Callable, Container, Sequence
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
    Check if a Python object is a ``bool``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, bool):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,bool.__name__,type(v).__name__))

# def check_is_bytearray(v: object, vname: str) -> None:
#     generic_check_isinstance(v, vname, bytearray)

# def check_is_bytes(v: object, vname: str) -> None:
#     generic_check_isinstance(v, vname, bytes)

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

# def check_is_complex(v: object, vname: str) -> None:
#     generic_check_isinstance(v, vname, complex)

def check_is_dict(v: object, vname: str) -> None:
    """
    Check if a Python object is a ``dict``. Raise error if not.

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
    Check if a Python object is a ``float``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, float):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,float.__name__,type(v).__name__))

# def check_is_frozenset(v: object, vname: str) -> None:
#     generic_check_isinstance(v, vname, frozenset)

def check_is_int(v: object, vname: str) -> None:
    """
    Check if a Python object is an ``int``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, int):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,int.__name__,type(v).__name__))

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

def check_is_list(v: object, vname: str) -> None:
    """
    Check if a Python object is a ``list``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, list):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,list.__name__,type(v).__name__))

# def check_is_memoryview(v: object, vname: str) -> None:
#     generic_check_isinstance(v, vname, memoryview)

def check_is_range(v: object, vname: str) -> None:
    """
    Check if a Python object is a ``range``. Raise error if not.

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
    Check if a Python object is a ``set``. Raise error if not.

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
    Check if a Python object is a ``str``. Raise error if not.

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
    Check if a Python object is a ``tuple``. Raise error if not.

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
    Check if a Python object is a ``type``. Raise error if not.

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
    Check if a Python object is a ``numbers.Complex``. Raise error if not.

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
    Check if a Python object is a ``numbers.Integral``. Raise error if not.

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
    Check if a Python object is a ``numbers.Number``. Raise error if not.

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
    Check if a Python object is a ``numbers.Real``. Raise error if not.

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
    Check if a Python object is a ``typing.Callable``. Raise error if not.

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
    Check if a Python object is a ``typing.Container``. Raise error if not.

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
    Check if a Python object is a ``typing.Sequence``. Raise error if not.

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
