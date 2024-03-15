"""
Module defining comma separated value I/O interfaces and associated error
checking routines.
"""

__all__ = [
    "CSVInputOutput",
    "check_is_CSVInputOutput",
]

from abc import ABCMeta, abstractmethod

class CSVInputOutput(
        metaclass = ABCMeta,
    ):
    """
    Abstract class for defining CSV input/output functionality.

    This abstract class defines two functions with the following purposes:

    - ``to_csv`` - write an object to a csv file.
    - ``from_csv`` - load an object from a csv file.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################

    ####################### File I/O #######################
    @abstractmethod
    def to_csv(
            self, 
            filename: str,
            **kwargs: dict
        ) -> None:
        """
        Write an object to a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name to which to write.
        kwargs : dict
            Additional keyword arguments to use for dictating export to a CSV.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ####################### File I/O #######################
    @classmethod
    @abstractmethod
    def from_csv(
            cls, 
            filename: str,
            **kwargs: dict
        ) -> 'CSVInputOutput':
        """
        Read an object from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : CSVInputOutput
            An object read from a CSV file.
        """
        raise NotImplementedError("class method is abstract")



################################## Utilities ###################################
def check_is_CSVInputOutput(v: object, vname: str) -> None:
    """
    Check if object is of type CSVInputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, CSVInputOutput):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,CSVInputOutput.__name__,type(v).__name__))
