"""
Module defining interfaces for selection problems
"""

# list of public objects in this module
__all__ = [
    "SelectionProblemType",
    "check_is_SelectionProblemType"
]

# imports
from abc import abstractmethod
from numbers import Integral
from typing import Callable, Tuple
import numpy
from pybrops.opt.prob.ProblemType import ProblemType

class SelectionProblemType(ProblemType):
    """
    Basal interface for all genotype selection problem specifications.

    All selection optimization problems have the form:

    .. math::

        \\min_{\\mathbf{x}} \\mathbf{w_{obj} \\odot T_{obj}(L(x))}

    Such that:

    .. math::

        \\mathbf{w_{ineqcv} \\odot T_{ineqcv}(L(x)) \\leq 0}

        \\mathbf{w_{eqcv} \\odot T_{eqcv}(L(x)) = 0}
    
    Where:

        - :math:`\\mathbf{x}` is a selection decision vector.
        - :math:`L(\\cdot)` is a latent vector encoding function transforming 
            the decision vector into a latent space. The ``latentfn`` function 
            defines this function in this interface.
        - :math:`w_{obj}` is an objective function weight vector.
        - :math:`T_{obj}(\\cdot)` is a function transforming a latent space 
            vector to an objective space vector.
        - :math:`w_{ineqcv}` is an inequality constraint violation function 
            weight vector.
        - :math:`T_{ineqcv}(\\cdot)` is a function transforming a latent space 
            vector to an inequality constraint violation vector.
        - :math:`w_{eqcv}` is an equality constraint violation function 
            weight vector.
        - :math:`T_{eqcv}(\\cdot)` is a function transforming a latent space 
            vector to an equality constraint violation vector.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SelectionProblemType.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(SelectionProblemType, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    @abstractmethod
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def obj_trans(self) -> Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to objective function values."""
        raise NotImplementedError("property is abstract")
    @obj_trans.setter
    @abstractmethod
    def obj_trans(self, value: Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]) -> None:
        """Set latent space to objective space transformation function."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def obj_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to objective space transformation function."""
        raise NotImplementedError("property is abstract")
    @obj_trans_kwargs.setter
    @abstractmethod
    def obj_trans_kwargs(self, value: dict) -> None:
        """Set keyword arguments for the latent space to objective space transformation function."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def ineqcv_trans(self) -> Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to inequality constraint violation values."""
        raise NotImplementedError("property is abstract")
    @ineqcv_trans.setter
    @abstractmethod
    def ineqcv_trans(self, value: Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]) -> None:
        """Set latent space to inequality constraint violation transformation function."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def ineqcv_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to inequality constraint violation transformation function."""
        raise NotImplementedError("property is abstract")
    @ineqcv_trans_kwargs.setter
    @abstractmethod
    def ineqcv_trans_kwargs(self, value: dict) -> None:
        """Set keyword arguments for the latent space to inequality constraint violation transformation function."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def eqcv_trans(self) -> Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to equality constraint violation values."""
        raise NotImplementedError("property is abstract")
    @eqcv_trans.setter
    @abstractmethod
    def eqcv_trans(self, value: Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]) -> None:
        """Set latent space to equality constraint violation transformation function."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def eqcv_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to equality constraint violation transformation function."""
        raise NotImplementedError("property is abstract")
    @eqcv_trans_kwargs.setter
    @abstractmethod
    def eqcv_trans_kwargs(self, value: dict) -> None:
        """Set keyword arguments for the latent space to equality constraint violation transformation function."""
        raise NotImplementedError("property is abstract")
    
    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    @abstractmethod
    def evalfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]:
        """
        Evaluate a candidate solution for the given Problem.
        
        This calculates three vectors:

        .. math::

            \\mathbf{v_{obj}} = \\mathbf{w_{obj} \\odot T_{obj}(L(x))} \\
            \\mathbf{v_{ineqcv}} = \\mathbf{w_{ineqcv} \\odot T_{ineqcv}(L(x))} \\
            \\mathbf{v_{eqcv}} = \\mathbf{w_{eqcv} \\odot T_{eqcv}(L(x))}

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
        out : tuple
            A tuple ``(obj, ineqcv, eqcv)``.
            
            Where:
            
            - ``obj`` is a numpy.ndarray of shape ``(nobj,)`` that contains 
                objective function evaluations.
                This is equivalent to :math:`\\mathbf{v_{obj}}`
            - ``ineqcv`` is a numpy.ndarray of shape ``(nineqcv,)`` that contains 
                inequality constraint violation values.
                This is equivalent to :math:`\\mathbf{v_{ineqcv}}`
            - ``eqcv`` is a numpy.ndarray of shape ``(neqcv,)`` that contains 
                equality constraint violation values.
                This is equivalent to :math:`\\mathbf{v_{eqcv}}`
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def latentfn(
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
def check_is_SelectionProblemType(v: object, vname: str) -> None:
    """
    Check if object is of type SelectionProblemType, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SelectionProblemType):
        raise TypeError("'{0}' must be of type SelectionProblemType.".format(vname))
