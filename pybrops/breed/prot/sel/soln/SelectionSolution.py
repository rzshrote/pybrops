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
