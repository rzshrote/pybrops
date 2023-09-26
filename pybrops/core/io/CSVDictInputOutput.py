"""
Module defining comma separated value I/O interfaces and associated error
checking routines.
"""

__all__ = [
    "CSVDictInputOutput",
    "check_is_CSVDictInputOutput",
]

from abc import ABCMeta, abstractmethod
from typing import Dict


class CSVDictInputOutput(metaclass=ABCMeta):
    """
    Abstract class for defining CSV input/output functionality.

    This abstract class defines two functions with the following purposes:

    - ``to_csv_dict`` - write an object to a csv file.
    - ``from_csv_dict`` - load an object from a csv file.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################

    ####################### File I/O #######################
    @abstractmethod
    def to_csv_dict(
            self, 
            filenames: Dict[str,str],
            **kwargs: dict
        ) -> None:
        """
        Write an object to a set of CSV files specified by values in a ``dict``.

        Parameters
        ----------
        filenames : dict
            Dictionary of CSV file names to which to write.
        kwargs : dict
            Additional keyword arguments to use for dictating export to a CSV.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ####################### File I/O #######################
    @classmethod
    @abstractmethod
    def from_csv_dict(
            cls, 
            filenames: Dict[str,str],
            **kwargs: dict
        ) -> 'CSVDictInputOutput':
        """
        Read an object from a set of CSV files specified by values in a ``dict``.

        Parameters
        ----------
        filenames : str
            Dictionary of CSV file names from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : CSVDictInputOutput
            An object read from a set of CSV files.
        """
        raise NotImplementedError("class method is abstract")



################################## Utilities ###################################
def check_is_CSVDictInputOutput(v: object, vname: str) -> None:
    """
    Check if object is of type CSVDictInputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, CSVDictInputOutput):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,CSVDictInputOutput.__name__,type(v).__name__))
