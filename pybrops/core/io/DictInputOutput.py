"""
Module defining ``dict`` I/O interfaces and associated error
checking routines.
"""

__all__ = [
    "DictInputOutput",
    "check_is_DictInputOutput",
]

from abc import ABCMeta, abstractmethod

class DictInputOutput(
        metaclass = ABCMeta,
    ):
    """
    Abstract class for defining ``dict`` input/output functionality.

    This abstract class defines two functions with the following purposes:

    - ``to_dict`` - export an object to a ``dict``.
    - ``from_dict`` - load an object from a ``dict``.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################

    ####################### File I/O #######################
    @abstractmethod
    def to_dict(
            self, 
            **kwargs: dict
        ) -> dict:
        """
        Export an object to a ``dict``.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            ``dict``.
        
        Returns
        -------
        out : dict
            An output Python dictionary.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ####################### File I/O #######################
    @classmethod
    @abstractmethod
    def from_dict(
            cls, 
            dic: dict,
            **kwargs: dict
        ) -> 'DictInputOutput':
        """
        Read an object from a ``dict``.

        Parameters
        ----------
        dic : dict
            Python dictionary from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            ``dict``.

        Returns
        -------
        out : DictInputOutput
            An object read from a ``dict``.
        """
        raise NotImplementedError("class method is abstract")



################################## Utilities ###################################
def check_is_DictInputOutput(v: object, vname: str) -> None:
    """
    Check if object is of type DictInputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DictInputOutput):
        raise TypeError(
            "variable '{0}' must be a of type '{1}' but received type '{2}'".format(
                vname,
                DictInputOutput.__name__,
                type(v).__name__
            )
        )
