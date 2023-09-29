"""
Module containing subroutines to check Pandas values.
"""

__all__ = [
    "check_DataFrame_has_column",
    "check_DataFrame_has_columns",
]
    
from numbers import Integral
import pandas


##################################################
########### DataFrame check functions ############
##################################################

def check_DataFrame_has_column(v: pandas.DataFrame, vname: str, column: str) -> None:
    """
    Check if a ``pandas.DataFrame`` contains a required column.

    Parameters
    ----------
    v : pandas.DataFrame
        Input ``pandas.DataFrame``.
    vname : str
        Name of the ``pandas.DataFrame`` variable.
    column : object
        A required column that the ``pandas.DataFrame`` must contain.
    """
    if column not in v.columns:
        raise ValueError("DataFrame '{0}' must have column '{1}'".format(vname,column))

def check_DataFrame_has_columns(v: pandas.DataFrame, vname: str, *args: tuple[str]) -> None:
    """
    Check if a ``pandas.DataFrame`` contains all required columns.

    Parameters
    ----------
    v : pandas.DataFrame
        Input ``pandas.DataFrame``.
    vname : str
        Name of the ``pandas.DataFrame`` variable.
    args : tuple
        Required columns that the ``pandas.DataFrame`` must contain.
    """
    for column in args:
        if column not in v.columns:
            raise ValueError("DataFrame '{0}' must have column '{1}'".format(vname,column))

def check_DataFrame_has_column_index(v: pandas.DataFrame, vname: str, ix: Integral) -> None:
    """
    Check if a ``pandas.DataFrame`` contains a required column index.

    Parameters
    ----------
    v : pandas.DataFrame
        Input ``pandas.DataFrame``.
    vname : str
        Name of the ``pandas.DataFrame`` variable.
    ix : object
        A required column index that the ``pandas.DataFrame`` must contain.
    """
    if (ix < 0) or (len(v.columns) <= ix):
        raise ValueError("DataFrame '{0}' must have column index '{1}'".format(vname,ix))

def check_DataFrame_has_column_indices(v: pandas.DataFrame, vname: str, *args: tuple[Integral]) -> None:
    """
    Check if a ``pandas.DataFrame`` contains a required column indices.

    Parameters
    ----------
    v : pandas.DataFrame
        Input ``pandas.DataFrame``.
    vname : str
        Name of the ``pandas.DataFrame`` variable.
    args : tuple
        A required column that the ``pandas.DataFrame`` must contain.
    """
    for ix in args:
        if (ix < 0) or (len(v.columns) <= ix):
            raise ValueError("DataFrame '{0}' must have column index '{1}'".format(vname,ix))

##################################################
############# Series check functions #############
##################################################

def check_Series_has_value(v: pandas.Series, vname: str, value: object) -> None:
    """
    Check if a ``pandas.Series`` contains a required value.

    Parameters
    ----------
    v : Series
        Input ``pandas.Series``.
    vname : str
        Name of the ``pandas.Series`` variable.
    value : object
        A required value that the ``pandas.Series`` must contain.
    """
    if value not in v:
        raise ValueError("pandas.Series '{0}' must have value '{1}'".format(vname,value))

def check_Series_has_values(v: pandas.Series, vname: str, *args: tuple[str]) -> None:
    """
    Check if a ``pandas.Series`` contains all required values.

    Parameters
    ----------
    v : Series
        Input ``pandas.Series``.
    vname : str
        Name of the ``pandas.Series`` variable.
    args : tuple
        Required values that the ``pandas.Series`` must contain.
    """
    for value in args:
        if value not in v:
            raise ValueError("pandas.Series '{0}' must have value '{1}'".format(vname,value))

def check_Series_has_index(v: pandas.Series, vname: str, ix: Integral) -> None:
    """
    Check if a ``pandas.Series`` contains a required index.

    Parameters
    ----------
    v : Series
        Input ``pandas.Series``.
    vname : str
        Name of the ``pandas.Series`` variable.
    ix : object
        A required index that the ``pandas.Series`` must contain.
    """
    if (ix < 0) or (len(v) <= ix):
        raise ValueError("pandas.Series '{0}' must have index '{1}'".format(vname,ix))

def check_Series_has_indices(v: pandas.Series, vname: str, *args: tuple[Integral]) -> None:
    """
    Check if a ``pandas.Series`` contains required indices.

    Parameters
    ----------
    v : Series
        Input ``pandas.Series``.
    vname : str
        Name of the ``pandas.Series`` variable.
    args : tuple
        Required indices that the ``pandas.Series`` must contain.
    """
    for ix in args:
        if (ix < 0) or (len(v) <= ix):
            raise ValueError("pandas.Series '{0}' must have index '{1}'".format(vname,ix))

