"""
Module containing error subroutines related to Python object attributes.
"""

__all__ = [
    "error_readonly",
    "check_is_callable",
    "check_is_iterable",
]

### read/write ###
def error_readonly(vname: str):
    """
    Subroutine to raise a read-only AttributeError with a custom error message.

    Parameters
    ----------
    vname : str
        Name of a variable.
    """
    raise AttributeError("variable '{0}' is read-only".format(vname))

################################################################################
########################### attribute check functions ##########################
################################################################################
def check_is_callable(v: object, vname: str) -> None:
    """
    Subroutine to check whether a Python object is callable.
    If the object is not callable, raise an AttributeError with a custom error
    message.

    Parameters
    ----------
    v : object
        Any Python object variable.
    vname : str
        Name of the Python object variable.
    """
    if not hasattr(v, "__call__"):
        raise AttributeError("variable '{0}' must be callable (have the '__call__' attribute)".format(vname))

def check_is_iterable(v: object, vname: str) -> None:
    """
    Subroutine to check whether a Python object is iterable.
    If the object is not iterable, raise an AttributeError with a custom error
    message.

    Parameters
    ----------
    v : object
        Any Python object variable.
    vname : str
        Name of the Python object variable.
    """
    if not hasattr(v, "__iter__"):
        raise AttributeError("variable '{0}' must be iterable (have the '__iter__' attribute)".format(vname))
