"""
Module containing selection solution interfaces
"""

__all__ = [
    "SelectionSolution"
]

from abc import ABCMeta
from pybrops.opt.soln.Solution import Solution


class SelectionSolution(Solution,metaclass=ABCMeta):
    """
    Class representing selection solutions.
    """
    # interface identical to Solution interface
    pass



################################## Utilities ###################################
def check_is_SelectionSolution(v: object, vname: str) -> None:
    """
    Check if object is of type SelectionSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SelectionSolution):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,SelectionSolution.__name__,type(v).__name__))
