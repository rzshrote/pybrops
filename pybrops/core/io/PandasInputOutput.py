"""
Module defining pandas.DataFrame I/O interfaces and associated error
checking routines.
"""

__all__ = [
    "PandasInputOutput",
    "check_is_PandasInputOutput",
]

from abc import ABCMeta, abstractmethod

import pandas


class PandasInputOutput(metaclass=ABCMeta):
    """
    Abstract class for defining pandas.DataFrame input/output functionality.

    This abstract class defines two functions with the following purposes:

    - ``to_pandas`` - export an object to a pandas.DataFrame.
    - ``from_pandas`` - load an object from a pandas.DataFrame.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################

    ####################### File I/O #######################
    @abstractmethod
    def to_pandas(
            self, 
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Export an object to a pandas.DataFrame.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            pandas.DataFrame.
        
        Returns
        -------
        out : pandas.DataFrame
            An output dataframe.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ####################### File I/O #######################
    @classmethod
    @abstractmethod
    def from_pandas(
            cls, 
            df: pandas.DataFrame,
            **kwargs: dict
        ) -> 'PandasInputOutput':
        """
        Read an object from a pandas.DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            pandas.DataFrame.

        Returns
        -------
        out : PandasInputOutput
            An object read from a pandas.DataFrame.
        """
        raise NotImplementedError("class method is abstract")



################################## Utilities ###################################
def check_is_PandasInputOutput(v: object, vname: str) -> None:
    """
    Check if object is of type PandasInputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, PandasInputOutput):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,PandasInputOutput.__name__,type(v).__name__))
