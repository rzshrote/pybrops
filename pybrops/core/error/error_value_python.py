from typing import Tuple

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
        raise ValueError("the length of '{0}' is not equal to {1}: received length {2}".format(vname,vlen,len(v)))

def check_len_eq(v: object, vname: str, w: object, wname: str) -> None:
    """
    Check is two objects have the same lengths.
    Raise error if object lengths are not equivalent.

    Parameters
    ----------
    v : object
        The first Python object for which to compare lengths.
    vname : str
        Name of the first Python object for which to compare lengths.
    w : object
        The second Python object for which to compare lengths.
    wname : str
        Name of the second Python object for which to compare lengths.
    """
    if len(v) != len(w):
        raise ValueError("the lengths of '{0}' and '{1}' are not equivalent: received lengths {2} and {3}, respectively".format(vname,wname,len(v),len(w)))

def check_all_equal(v: object, vname: str) -> None:
    """
    Check whether an object's elements are all equivalent.

    Parameters
    ----------
    v : object
        A Python object.
    vname : object
        Name of the Python object.
    """
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
        A Python object.
    vname : str
        Name of the Python object for use in the error message.
    """
    if v < 0:
        raise ValueError("variable '{0}' must be positive".format(vname))

def check_is_eq(v: object, vname: str, value: object) -> None:
    """
    Check if a Python object is equivalent to another Python object.

    Parameters
    ----------
    v : object
        A Python object.
    vname : str
        Name of the Python object for use in the error message.
    value : object
        Expected value of the input Python object.
    """
    if v != value:
        raise ValueError("variable '{0}' is not equal to {1}".format(vname, value))

def check_is_gt(v: object, vname: str, value: object) -> None:
    """
    Check if a Python object is greater than another Python object.

    Parameters
    ----------
    v : object
        A Python object.
    vname : str
        Name of the Python object for use in the error message.
    value : object
        Lower value of the input Python object.
    """
    if v <= value:
        raise ValueError("variable '{0}' is not greater than {1}".format(vname, value))

def check_is_gteq(v: object, vname: str, value: object) -> None:
    """
    Check if a Python object is greater than or equal to another Python object.

    Parameters
    ----------
    v : object
        A Python object.
    vname : str
        Name of the Python object for use in the error message.
    value : object
        Lower value of the input Python object.
    """
    if v < value:
        raise ValueError("variable '{0}' is not greater than or equal to {1}".format(vname, value))

def check_is_lt(v: object, vname: str, value: object) -> None:
    """
    Check if a Python object is less than another Python object.

    Parameters
    ----------
    v : object
        A Python object.
    vname : str
        Name of the Python object for use in the error message.
    value : object
        Upper value of the input Python object.
    """
    if v >= value:
        raise ValueError("variable '{0}' is not less than {1}".format(vname, value))

def check_is_lteq(v: object, vname: str, value: object) -> None:
    """
    Check if a Python object is less than or equal to another Python object.

    Parameters
    ----------
    v : object
        A Python object.
    vname : str
        Name of the Python object for use in the error message.
    value : object
        Upper value of the input Python object.
    """
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
def check_keys_in_dict(v: dict, vname: str, *args: Tuple[object,...]) -> None:
    """
    Check if a set of keys can be found among the keys of a provided ``dict``.
    If any of the keys cannot be found within the dict, then raise an error.

    Parameters
    ----------
    v : dict
        Dictionary to check for keys.
    vname : str
        Name assigned to the dictionary.
    args : Tuple[object,...]
        A tuple of keys for which to search in the input dictionary.
    """
    keys = v.keys()
    if any(e not in keys for e in args):
        raise ValueError("dict '{0}' must have keys: {1}".format(vname, args))

def check_keys_in_dict_all_type(v: dict, vname: str, vtype: type) -> None:
    """
    Check if all keys in a dictionary are of a specific type.

    Parameters
    ----------
    v : dict
        Dictionary in which to check keys.
    vname : str
        Name assigned to the dictionary.
    vtype : type
        Expected type of the keys in the dictionary.
    """
    if any(not isinstance(e, vtype) for e in v.keys()):
        raise ValueError("not all keys in dict '{0}' are of type {1}".format(vname, str(vtype)))

def check_values_in_dict_all_type(v: dict, vname: str, vtype: type) -> None:
    """
    Check if all values in a dictionary have a specific type.

    Parameters
    ----------
    v : dict
        Dictionary in which to check values.
    vname : str
        Name assigned to the dictionary.
    vtype : type
        Expected type of the values in the dictionary.
    """
    if any(not isinstance(e, vtype) for e in v.values()):
        raise ValueError("not all values in dict '{0}' are of type {1}".format(vname, str(vtype)))

def check_values_in_dict_equal_len(v: dict, vname: str) -> None:
    """
    Check if all values in a ``dict`` have a equal lengths.

    Parameters
    ----------
    v : dict
        Dictionary in which to check values.
    vname : str
        Name assigned to the dictionary.
    """
    viter = iter(v.values())
    try:
        e0 = next(viter)
    except StopIteration:
        return
    l0 = len(e0)
    if any(len(e) != l0 for e in viter):
        raise ValueError("not all values in dict '{0}' have equal length == {1}".format(vname, l0))

def check_values_in_dict_len(v: dict, vname: str, vlen: int) -> None:
    """
    Check if all values in a ``dict`` have a specific length.

    Parameters
    ----------
    v : dict
        Dictionary in which to check values.
    vname : str
        Name assigned to the dictionary.
    vlen : int
        Expected length for each value in the dictionary.
    """
    viter = iter(v.values())
    try:
        e0 = next(viter)
    except StopIteration:
        raise ValueError("dict '{0}' is empty".format(vname))
    if any(len(e) != vlen for e in viter):
        raise ValueError("not all values in dict '{0}' have length == {1}".format(vname, vlen))
