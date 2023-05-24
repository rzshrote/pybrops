import numbers
from typing import Any

################################################################################
############################### check functions ################################
################################################################################
def check_is_not_None(v: object, vname: str) -> None:
    """
    Check if an object is not None.
    Raise error if object is None.

    Parameters
    ----------
    v : object
        Any Python object
    vname : str
        Name of the Python object for use in the error message.
    """
    if v is None:
        raise ValueError("variable '{0}' is not 'None'".format(vname))

def check_len(v: object, vname: str, vlen: int) -> None:
    """
    Check if an object has a length equal to a provided value.
    Raise error if object length is not equal to provided value.

    Parameters
    ----------
    v : object
        Any Python object
    vname : str
        Name of the Python object for use in the error message.
    vlen : int
        Expected length of the Python object.
    """
    if len(v) != vlen:
        raise ValueError("the length of '{0}' is not equal to {1}".format(vname,vlen))

def check_len_eq(v: object, vname: str, w: Any, wname: str) -> None:
    if len(v) != len(w):
        raise ValueError("the lengths of '{0}' and '{1}' are not equivalent".format(vname, wname))

def check_all_equal(v: object, vname: str) -> None:
    viter = iter(v)
    try:
        e0 = next(viter)
    except StopIteration:
        return # length is 1, therefore all elements are equivalent
    if any(e0 != e for e in viter):
        raise ValueError("not all elements in {0} are equal to {1}".format(vname, e0))

def check_is_positive(v: object, vname: str) -> None:
    """
    Check if a Python object is positive. Raise error if it is not.

    Parameters
    ----------
    v : object
        Any Python object,
    vname : str
        Name of the Python object for use in the error message.
    """
    if v < 0:
        raise ValueError("variable '{0}' must be positive".format(vname))

def check_is_eq(v: object, vname: str, value: object) -> None:
    """Raise error if ``obj`` is not equal to ``value``."""
    if v != value:
        raise ValueError("variable '{0}' is not equal to {1}".format(vname, value))

def check_is_gt(v: object, vname: str, value: object):
    """Raise error if ``obj`` is not greater than ``value``."""
    if v <= value:
        raise ValueError("variable '{0}' is not greater than {1}".format(vname, value))

def check_is_gteq(v: object, vname: str, value: object):
    """Raise error if ``obj`` is not greater than ``value``."""
    if v < value:
        raise ValueError("variable '{0}' is not greater than {1}".format(vname, value))

def check_is_lt(v: object, vname: str, value: object):
    """Raise error if ``v`` is not less than ``value``."""
    if v >= value:
        raise ValueError("variable '{0}' is not less than {1}".format(vname, value))

def check_is_lteq(v: object, vname: str, value: object):
    """Raise error if ``v`` is not less than ``value``."""
    if v > value:
        raise ValueError("variable '{0}' is not less than or equal to {1}".format(vname, value))

def check_is_in_interval(v: object, vname: str, vmin: object, vmax: object) -> None:
    """
    Check if a value is in a specified interval.
    
    Parameters
    ----------
    v : object
        Input numerical value.
    vname : str
        Name of the variable of which to test.
    vmin : object
        Minimum value for the input value.
    vmax : object
        Maximum value for the input value.
    
    Returns
    -------
    out : outtype
        outdesc
    """
    if (v < vmin) or (v > vmax):
        raise ValueError("variable '{0}' is not in interval [{1}, {2}]".format(vname, vmin, vmax))

def check_str_value(v: str, vname: str, *args: tuple) -> None:
    """
    Check if a string has an accepted value:

    Parameters
    ----------
    v : str
        Input string.
    vname : str
        Name of the variable for the string.
    args : tuple
        Acceptable values for the string.
    """
    if v not in args:
        raise ValueError("string '{0}' must be one of: {1}".format(vname, args))

##################################################
########### Dictionary check functions ###########
##################################################
def check_keys_in_dict(v, vname, *args):
    keys = v.keys()
    if any(e not in keys for e in args):
        raise ValueError("dict '{0}' must have keys: {1}".format(vname, args))

def check_keys_in_dict_all_type(v, vname, vtype):
    if any(not isinstance(e, vtype) for e in v.keys()):
        raise ValueError("not all keys in dict '{0}' are of type {1}".format(vname, str(vtype)))

def check_values_in_dict_all_type(v, vname, vtype):
    if any(not isinstance(e, vtype) for e in v.values()):
        raise ValueError("not all values in dict '{0}' are of type {1}".format(vname, str(vtype)))

def check_values_in_dict_equal_len(v: object, vname: str) -> None:
    viter = iter(v.values())
    try:
        e0 = next(viter)
    except StopIteration:
        return
    l0 = len(e0)
    if any(len(e) != l0 for e in viter):
        raise ValueError("not all values in dict '{0}' have equal length == {1}".format(vname, l0))

def check_values_in_dict_len(v, vname, l):
    viter = iter(v.values())
    try:
        e0 = next(viter)
    except StopIteration:
        raise ValueError("dict '{0}' is empty".format(vname))
    if any(len(e) != l for e in viter):
        raise ValueError("not all values in dict '{0}' have length == {1}".format(vname, l))
