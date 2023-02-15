import math
import numbers
from typing import Any
from typing import Union
import inspect

################################################################################
###################### basic inheritance check functions #######################
################################################################################
def check_inherits(v, vname, vtype):
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
def check_isinstance(v: Any, vname: str, vtype: Union[type,tuple]) -> None:
    """
    Generic check type function.

    Parameters
    ----------
    v : Any
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

def check_is_bool(v: Any, vname: str) -> None:
    """
    Check if a Python object is a bool. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, bool):
        raise TypeError("variable '{0}' must be of type 'bool'".format(vname))

# def check_is_bytearray(v: Any, vname: str) -> None:
#     generic_check_isinstance(v, vname, bytearray)

# def check_is_bytes(v: Any, vname: str) -> None:
#     generic_check_isinstance(v, vname, bytes)

def check_is_class(v: Any, vname: str) -> None:
    if not inspect.isclass(v):
        raise TypeError("variable '{0}' must be a class name".format(vname))

# def check_is_complex(v: Any, vname: str) -> None:
#     generic_check_isinstance(v, vname, complex)

def check_is_dict(v: Any, vname: str) -> None:
    """
    Check if a Python object is a dict. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, dict):
        raise TypeError("variable '{0}' must be of type 'dict'".format(vname))

def check_is_float(v: Any, vname: str) -> None:
    """
    Check if a Python object is a float. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, float):
        raise TypeError("variable '{0}' must be of type 'float'".format(vname))

# def check_is_frozenset(v: Any, vname: str) -> None:
#     generic_check_isinstance(v, vname, frozenset)

def check_is_int(v: Any, vname: str) -> None:
    """
    Check if a Python object is a int. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, int):
        raise TypeError("variable '{0}' must be of type 'int'".format(vname))

def check_is_int_or_inf(v: Any, vname: str) -> None:
    """
    Check if a Python object is an int or if it is infinite. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if (v != math.inf) and (not isinstance(v, int)):
        raise TypeError("variable '{0}' must be of type 'int' or be infinite".format(vname))

def check_is_list(v: Any, vname: str) -> None:
    """
    Check if a Python object is a list. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, list):
        raise TypeError("variable '{0}' must be of type 'list'".format(vname))

# def check_is_memoryview(v: Any, vname: str) -> None:
#     generic_check_isinstance(v, vname, memoryview)

def check_is_range(v: Any, vname: str) -> None:
    """
    Check if a Python object is a range. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, range):
        raise TypeError("variable '{0}' must be of type 'range'".format(vname))

def check_is_set(v: Any, vname: str) -> None:
    """
    Check if a Python object is a set. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, set):
        raise TypeError("variable '{0}' must be of type 'set'".format(vname))

def check_is_str(v: Any, vname: str) -> None:
    """
    Check if a Python object is a str. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, str):
        raise TypeError("variable '{0}' must be of type 'str'".format(vname))

def check_is_tuple(v: Any, vname: str) -> None:
    """
    Check if a Python object is a tuple. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, tuple):
        raise TypeError("variable '{0}' must be of type 'tuple'".format(vname))

def check_is_type(v: Any, vname: str) -> None:
    """
    Check if a Python object is a type. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, type):
        raise TypeError("variable '{0}' must be of type 'type'".format(vname))

def check_is_Number(v: Any, vname: str) -> None:
    """
    Check if a Python object is a numbers.Number. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, numbers.Number):
        raise TypeError("variable '{0}' must be of type 'numbers.Number'".format(vname))

def check_is_Integral(v: Any, vname: str) -> None:
    """
    Check if a Python object is a numbers.Integral. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, numbers.Integral):
        raise TypeError("variable '{0}' must be of type 'numbers.Integral'".format(vname))

################################################################################
################ compound check functions for basic data types #################
################################################################################
def check_is_array_like(v: Any, vname: str) -> None:
    alattr = ("__len__","__iter__","__getitem__")
    for a in alattr:
        if not hasattr(v, a):
            raise TypeError("'{0}' must have attribute '{1}' to be array_like".format(vname,a))

def check_is_str_or_iterable(v: Any, vname: str) -> None:
    if not (isinstance(v, str) or hasattr(v, "__iter__")):
        raise TypeError("'{0}' must be of type str or have attribute __iter__.".format(vname))

def check_is_list_or_tuple(v: Any, vname: str) -> None:
    """
    Check if a Python object is a list or tuple. Raise error if not.

    Parameters
    ----------
    v : Any
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, (list,tuple)):
        raise TypeError("variable '{0}' must be of type 'list' or 'tuple'".format(vname))
