"""
Module defining DataFrame interfaces and associated error checking routines.
"""

from typing import Any


class DataFrame:
    """
    An abstract class for data frame wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Data frame row and column access
        2) Data extraction from a data frame
        3) Data frame conversion
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class DataFrame.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(DataFrame, self).__init__()

    def __copy__(self):
        """
        Make a shallow copy of the the dataframe.

        Returns
        -------
        out : DataFrame
        """
        raise NotImplementedError("method is abstract")

    def __deepcopy__(self, memo):
        """
        Make a deep copy of the dataframe.

        Parameters
        ----------
        memo : dict

        Returns
        -------
        out : DataFrame
        """
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def data(self) -> Any:
        """Access to raw data frame object."""
        raise NotImplementedError("property is abstract")
    @data.setter
    def data(self, value: Any) -> None:
        """Set dataframe"""
        raise NotImplementedError("property is abstract")
    @data.deleter
    def data(self) -> None:
        """Delete dataframe"""
        raise NotImplementedError("property is abstract")

    ################## Column attributes ###################
    @property
    def ncol(self) -> Any:
        """Number of columns."""
        raise NotImplementedError("property is abstract")
    @ncol.setter
    def ncol(self, value: Any) -> None:
        """Set number of columns"""
        raise NotImplementedError("property is abstract")
    @ncol.deleter
    def ncol(self) -> None:
        """Delete number of columns"""
        raise NotImplementedError("property is abstract")

    @property
    def col_axis(self) -> int:
        """Column axis index."""
        raise NotImplementedError("property is abstract")
    @col_axis.setter
    def col_axis(self, value: int) -> None:
        """Set column axis index"""
        raise NotImplementedError("property is abstract")
    @col_axis.deleter
    def col_axis(self) -> None:
        """Delete column axis index"""
        raise NotImplementedError("property is abstract")

    @property
    def col_dtype(self) -> Any:
        """Column data types."""
        raise NotImplementedError("property is abstract")
    @col_dtype.setter
    def col_dtype(self, value: Any) -> None:
        """Set column data types"""
        raise NotImplementedError("property is abstract")
    @col_dtype.deleter
    def col_dtype(self) -> None:
        """Delete column data types"""
        raise NotImplementedError("property is abstract")
    
    @property
    def col_name(self) -> Any:
        """Column names."""
        raise NotImplementedError("property is abstract")
    @col_name.setter
    def col_name(self, value: Any) -> None:
        """Set column names"""
        raise NotImplementedError("property is abstract")
    @col_name.deleter
    def col_name(self) -> None:
        """Delete column names"""
        raise NotImplementedError("property is abstract")
    
    @property
    def col_grp(self) -> Any:
        """Column groups used for classifying variables."""
        raise NotImplementedError("property is abstract")
    @col_grp.setter
    def col_grp(self, value: Any) -> None:
        """Set column groups"""
        raise NotImplementedError("property is abstract")
    @col_grp.deleter
    def col_grp(self) -> None:
        """Delete column groups"""
        raise NotImplementedError("property is abstract")

    #################### Row attributes ####################
    @property
    def nrow(self) -> int:
        """Number of rows."""
        raise NotImplementedError("property is abstract")
    @nrow.setter
    def nrow(self, value: int) -> None:
        """Set number of rows"""
        raise NotImplementedError("property is abstract")
    @nrow.deleter
    def nrow(self) -> None:
        """Delete number of rows"""
        raise NotImplementedError("property is abstract")

    @property
    def row_axis(self) -> int:
        """Row axis index."""
        raise NotImplementedError("property is abstract")
    @row_axis.setter
    def row_axis(self, value: int) -> None:
        """Set row axis index"""
        raise NotImplementedError("property is abstract")
    @row_axis.deleter
    def row_axis(self) -> None:
        """Delete row axis index"""
        raise NotImplementedError("property is abstract")

    @property
    def row_name(self) -> Any:
        """Row names."""
        raise NotImplementedError("property is abstract")
    @row_name.setter
    def row_name(self, value: Any) -> None:
        """Set row names"""
        raise NotImplementedError("property is abstract")
    @row_name.deleter
    def row_name(self) -> None:
        """Delete row names"""
        raise NotImplementedError("property is abstract")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def col_data(self, **kwargs: dict):
        """
        Get a column's data from the dataframe.
        """
        raise NotImplementedError("method is abstract")

    def to_pandas_df(self, **kwargs: dict):
        """
        Get dataframe as a pandas.DataFrame.

        Returns
        -------
        out : pandas.DataFrame
            DataFrame as a pandas.DataFrame.
        """
        raise NotImplementedError("method is abstract")

    def to_dict(self, **kwargs: dict):
        """
        Get dataframe as a dictionary of numpy.ndarray's.

        Returns
        -------
        out : dict
            DataFrame as a dictionary of numpy.ndarray's.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_DataFrame(v: Any) -> bool:
    """
    Determine whether an object is a DataFrame.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DataFrame object instance.
    """
    return isinstance(v, DataFrame)

def check_is_DataFrame(v: Any, vname: str) -> None:
    """
    Check if object is of type DataFrame. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DataFrame):
        raise TypeError("variable '{0}' must be a DataFrame".format(vname))
