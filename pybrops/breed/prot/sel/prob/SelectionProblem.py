"""
Module defining interfaces for selection problems
"""

# list of public objects in this module
__all__ = [
    "SelectionProblem",
    "check_is_SelectionProblem"
]

# imports
from typing import Callable, Tuple
import numpy
from pybrops.opt.prob.Problem import Problem

class SelectionProblem(Problem):
    """
    docstring for SelectionProblem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(SelectionProblem, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def encode_trans(self) -> Callable[[numpy.ndarray,dict],Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]]:
        """Function which transforms outputs from ``encodefn`` to a tuple ``(obj,ineqcv,eqcv)``."""
        raise NotImplementedError("property is abstract")
    @encode_trans.setter
    def encode_trans(self, value: Callable[[numpy.ndarray,dict],Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]]) -> None:
        """Set ``encodefn`` output transformation function."""
        raise NotImplementedError("property is abstract")
    @encode_trans.deleter
    def encode_trans(self) -> None:
        """Delete ``encodefn`` output transformation function."""
        raise NotImplementedError("property is abstract")
    
    @property
    def encode_trans_kwargs(self) -> dict:
        """``encodefn`` output transformation function keyword arguments."""
        raise NotImplementedError("property is abstract")
    @encode_trans_kwargs.setter
    def encode_trans_kwargs(self, value: dict) -> None:
        """Set ``encodefn`` output transformation function keyword arguments."""
        raise NotImplementedError("property is abstract")
    @encode_trans_kwargs.deleter
    def encode_trans_kwargs(self) -> None:
        """Delete ``encodefn`` output transformation function keyword arguments."""
        raise NotImplementedError("property is abstract")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def encodefn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Encode a candidate solution for the given Problem into an ``l`` 
        dimensional latent evaluation space.
        
        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,)``.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(l,)`` containing latent evaluation values.
        """
        raise NotImplementedError("method is abstract")




################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SelectionProblem(v: object, vname: str) -> None:
    """
    Check if object is of type SelectionProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SelectionProblem):
        raise TypeError("'{0}' must be of type SelectionProblem.".format(vname))
