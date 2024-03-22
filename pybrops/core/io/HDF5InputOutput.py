"""
Module defining HDF5 I/O interfaces and associated error checking routines.
"""

__all__ = [
    "HDF5InputOutput",
    "check_is_HDF5InputOutput",
]

from abc import ABCMeta
from abc import abstractmethod
from typing import Optional
from typing import Union
import h5py
from pathlib import Path

class HDF5InputOutput(
        metaclass = ABCMeta
    ):
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
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str],
            overwrite: bool,
        ) -> None:
        """
        Write an object to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name to which to write.
            If ``Path``, an HDF5 file path to which to write.
            If ``h5py.File``, an opened HDF5 file to which to write.
        
        groupname : str, None
            If ``str``, an HDF5 group name under which object data is stored.
            If ``None``, object is written to the base HDF5 group.
        
        overwrite : bool
            Whether to overwrite values in an HDF5 file if a field already exists.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    @abstractmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str]
        ) -> 'HDF5InputOutput':
        """
        Read an object from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name from which to read.
            If ``Path``, an HDF5 file name from which to read.
            If ``h5py.File``, an opened HDF5 file from which to read. 
        
        groupname : str, None
            If ``str``, an HDF5 group name under which object data is stored.
            If ``None``, object is read from base HDF5 group.

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
