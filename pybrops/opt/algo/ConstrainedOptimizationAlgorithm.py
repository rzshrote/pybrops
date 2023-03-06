"""
Module defining interfaces and associated error checking routines for
constrained optimization algorithms.
"""

from numbers import Integral, Real
from typing import Any, Callable, Optional, Union
import numpy

class ConstrainedOptimizationAlgorithm:
    """
    An abstract class for optimization algorithms.

    The purpose of this abstract class is to provide functionality for:
        1) Optimization of objective functions.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedOptimizationAlgorithm.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(ConstrainedOptimizationAlgorithm, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def optimize(
            self, 
            evalfn: Callable, 
            ndecn: Integral, 
            decnspace: numpy.ndarray, 
            nobj: Integral,
            objwt: Union[Real,numpy.ndarray],
            nineqcv: Integral,
            ineqcvwt: Union[Real,numpy.ndarray],
            neqcv: Integral,
            eqcvwt: Union[Real,numpy.ndarray], 
            miscout: Optional[dict],
            **kwargs: dict
        ) -> tuple:
        """
        Optimize an objective function.

        Parameters
        ----------
        evalfn : Callable
            Evaluation function which to optimize. An evaluation function
            accepts a vector of decision variables and returns a tuple of 
            length 3 ``(obj, ineqcv, eqcv)`` containing objective function 
            value(s) (``obj``), inequality constraint violation value(s) 
            (``ineqcv``), and equality constraint violation value(s) 
            (``eqcv``).
        ndecn : int
            Number of decision variables in the search space.
            A vector is formed as decnspace**ndecn
        decnspace : numpy.ndarray
            Decision space in which the constrained optimization algorithm 
            searches.
        objwt : Real, numpy.ndarray
            Weight(s) applied to output(s) from the objective function.
            Weights are applied to the values held in ``obj``.
        ineqcvwt : Real, numpy.ndarray
            Weight(s) applied to output(s) from the inequality constraint 
            violation function. Weights are applied to the values held in 
            ``ineqcv``.
        eqcvwt : Real, numpy.ndarray
            Weight(s) applied to output(s) from the equality constraint 
            violation function. Weights are applied to the values held in 
            ``eqcv``.
        miscout : dict
            Miscellaneous output from the constrained optimizaiont algorithm.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : tuple
            A tuple of length 2 ``(soln, decn)`` containing the solution(s)
            (``soln``) and the corresponding decision variables for the 
            solution(s) (``decn``).
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_ConstrainedOptimizationAlgorithm(v: Any) -> bool:
    """
    Determine whether an object is a ConstrainedOptimizationAlgorithm.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a ConstrainedOptimizationAlgorithm object instance.
    """
    return isinstance(v, ConstrainedOptimizationAlgorithm)

def check_is_ConstrainedOptimizationAlgorithm(v: Any, vname: str) -> None:
    """
    Check if object is of type ConstrainedOptimizationAlgorithm. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, ConstrainedOptimizationAlgorithm):
        raise TypeError("variable '{0}' must be a ConstrainedOptimizationAlgorithm".format(vname))
