"""
Module containing error subroutines related to Python object attributes.
"""

from . import generic_default_cond

### read/write ###
def error_readonly(vname):
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
def check_is_callable(v, vname):
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
        raise AttributeError(
            "variable '{0}' must be callable (have the '__call__' attribute)".format(vname)
        )

def check_is_iterable(v, vname):
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
        raise AttributeError(
            "variable '{0}' must be callable (have the '__iter__' attribute)".format(vname)
        )

################################################################################
#################### conditional attribute check functions #####################
################################################################################
def cond_check_is_callable(v, vname, cond = generic_default_cond):
    """
    Subroutine to conditionally check whether a Python object is callable.
    If the object is not callable, raise an AttributeError with a custom error
    message.

    Parameters
    ----------
    v : object
        Any Python object variable.
    vname : str
        Name of the Python object variable.
    cond : function
        A Python function accepting ``v`` as an argument and returning
        ``True`` or ``False`` of whether to check if ``v`` is callable.
    """
    if cond(v):
        check_is_callable(v, vname)

def cond_check_is_iterable(v, vname, cond = generic_default_cond):
    """
    Subroutine to conditionally check whether a Python object is iterable.
    If the object is not iterable, raise an AttributeError with a custom error
    message.

    Parameters
    ----------
    v : object
        Any Python object variable.
    vname : str
        Name of the Python object variable.
    cond : function
        A Python function accepting ``v`` as an argument and returning
        ``True`` or ``False`` of whether to check if ``v`` is callable.
    """
    if cond(v):
        check_is_iterable(v, vname)
