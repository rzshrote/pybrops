"""
Module defining interfaces and associated error checking routines for
optimization algorithms.
"""

from typing import Callable

import numpy


class UnconstrainedOptimizationAlgorithm:
    """
    An abstract class for optimization algorithms.

    The purpose of this abstract class is to provide functionality for:
        1) Optimization of objective functions.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class OptimizationAlgorithm.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(UnconstrainedOptimizationAlgorithm, self).__init__()

    ############################## Object Methods ##############################
    def optimize(
            self, 
            objfn: Callable, 
            k: int, 
            sspace: numpy.ndarray, 
            objfn_wt : numpy.ndarray, 
            **kwargs: dict
        ):
        """
        Optimize an objective function.

        Parameters
        ----------
        objfn : callable
            Objective function which to optimize.
        k : int
            Number of decision variables in the search space.
            A vector is formed as sspace^k
        sspace : numpy.ndarray
            Search space that the OptimizationAlgorithm searches in.
        objfn_wt : numpy.ndarray
            Weight(s) applied to output(s) from the objfn.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : tuple
            A tuple of length 3 (soln, decn, misc)
        """
        raise NotImplementedError("method is abstract")

    # def optimize_vec(fn):
    #     raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_OptimizationAlgorithm(v: object, vname: str) -> None:
    """
    Check if object is of type OptimizationAlgorithm. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, UnconstrainedOptimizationAlgorithm):
        raise TypeError("variable '{0}' must be a OptimizationAlgorithm".format(vname))
