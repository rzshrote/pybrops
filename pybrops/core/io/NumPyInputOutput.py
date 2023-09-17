"""
Module defining numpy.ndarray I/O interfaces and associated error
checking routines.
"""

__all__ = [
    "NumPyInputOutput",
    "check_is_NumPyInputOutput",
]

from abc import ABCMeta, abstractmethod

import numpy

class NumPyInputOutput(metaclass=ABCMeta):
    """
    Abstract class for defining numpy.ndarray input/output functionality.

    This abstract class defines two functions with the following purposes:

    - ``to_numpy`` - export an object to a numpy.ndarray.
    - ``from_numpy`` - load an object from a numpy.ndarray.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################

    ####################### File I/O #######################
    @abstractmethod
    def to_numpy(
            self, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Export an object to a numpy.ndarray.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            numpy.ndarray.
        
        Returns
        -------
        out : numpy.ndarray
            An output dataframe.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ####################### File I/O #######################
    @classmethod
    @abstractmethod
    def from_numpy(
            cls, 
            arr: numpy.ndarray,
            **kwargs: dict
        ) -> 'NumPyInputOutput':
        """
        Read an object from a numpy.ndarray.

        Parameters
        ----------
        arr : numpy.ndarray
            NumPy array from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            numpy.ndarray.

        Returns
        -------
        out : NumPyInputOutput
            An object read from a numpy.ndarray.
        """
        raise NotImplementedError("class method is abstract")



################################## Utilities ###################################
def check_is_NumPyInputOutput(v: object, vname: str) -> None:
    """
    Check if object is of type NumPyInputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, NumPyInputOutput):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,NumPyInputOutput.__name__,type(v).__name__))
