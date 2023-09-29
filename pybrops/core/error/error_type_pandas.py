"""
Module containing subroutines to check Pandas types.
"""

__all__ = [
    "check_is_pandas_DataFrame",
]

from typing import Union
from pandas import DataFrame
import pandas

################################################################################
########################## isinstance check functions ##########################
################################################################################
def check_is_pandas_DataFrame(v: object, vname: str) -> None:
    """
    Check whether a Python object is a ``pandas.DataFrame``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, DataFrame):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,DataFrame.__name__,type(v).__name__))

def check_Series_all_type(v: pandas.Series, vname: str, vtype: Union[type,tuple]) -> None:
    if not all(isinstance(e, vtype) for e in v):
        raise TypeError(
            "pandas.Series '{0}' must have values all of type '{1}'".format(
                vname,
                vtype.__name__
            )
        )

