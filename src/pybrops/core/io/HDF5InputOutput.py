"""
Module defining HDF5 I/O interfaces and associated error checking routines.
"""

__all__ = [
    "HDF5InputOutput",
    "check_is_HDF5InputOutput",
]

from abc import ABCMeta, abstractmethod
from typing import Optional


class HDF5InputOutput(metaclass=ABCMeta):
    """
    Abstract class for defining HDF5 input/output functionality.

    This abstract class defines two functions with the following purposes:

    - ``to_hdf5`` - write an object to an HDF5 file.
    - ``from_hdf5`` - load an object from an HDF5 file.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################

    ################### Matrix File I/O ####################
    @abstractmethod
    def to_hdf5(
            self, 
            filename: str, 
            groupname: Optional[str]
        ) -> None:
        """
        Write an object to an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name to which to write.
        groupname : str or None
            HDF5 group name under which object data is stored.
            If None, object is written to the base HDF5 group.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    @abstractmethod
    def from_hdf5(
            cls, 
            filename: str, 
            groupname: Optional[str]
        ) -> 'HDF5InputOutput':
        """
        Read an object from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which object data is stored.
            If None, object is read from base HDF5 group.

        Returns
        -------
        out : HDF5InputOutput
            An object read from an HDF5 file.
        """
        raise NotImplementedError("class method is abstract")



################################## Utilities ###################################
def check_is_HDF5InputOutput(v: object, vname: str) -> None:
    """
    Check if object is of type HDF5InputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, HDF5InputOutput):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,HDF5InputOutput.__name__,type(v).__name__))
