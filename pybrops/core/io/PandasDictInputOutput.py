"""
Module defining pandas dictionary I/O interfaces and associated error
checking routines.
"""

__all__ = [
    "PandasDictInputOutput",
    "check_is_PandasDictInputOutput",
]

from abc import ABCMeta
from abc import abstractmethod
from typing import Dict

import pandas

class PandasDictInputOutput(
        metaclass = ABCMeta,
    ):
    """
    Abstract class for defining pandas dictionary input/output functionality.

    This abstract class defines two functions with the following purposes:

    - ``to_pandas_dict`` - export an object to a ``dict`` of ``pandas.DataFrame``.
    - ``from_pandas_dict`` - load an object from a ``dict`` of ``pandas.DataFrame``.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################

    ####################### File I/O #######################
    @abstractmethod
    def to_pandas_dict(
            self, 
            **kwargs: dict
        ) -> Dict[str,pandas.DataFrame]:
        """
        Export an object to a ``dict`` of ``pandas.DataFrame``.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            ``dict`` of ``pandas.DataFrame``.
        
        Returns
        -------
        out : dict
            An output dataframe.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ####################### File I/O #######################
    @classmethod
    @abstractmethod
    def from_pandas_dict(
            cls, 
            dic: Dict[str,pandas.DataFrame],
            **kwargs: dict
        ) -> 'PandasDictInputOutput':
        """
        Read an object from a ``dict`` of ``pandas.DataFrame``.

        Parameters
        ----------
        dic : dict
            Python dictionary containing ``pandas.DataFrame`` from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            ``dict`` of ``pandas.DataFrame``.

        Returns
        -------
        out : PandasDictInputOutput
            An object read from a ``dict`` of ``pandas.DataFrame``.
        """
        raise NotImplementedError("class method is abstract")



################################## Utilities ###################################
def check_is_PandasDictInputOutput(v: object, vname: str) -> None:
    """
    Check if object is of type ``PandasDictInputOutput``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, PandasDictInputOutput):
        raise TypeError(
            "variable '{0}' must be a of type '{1}' but received type '{2}'".format(
                vname,
                PandasDictInputOutput.__name__,
                type(v).__name__
            )
        )
