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
    def __init__(self, **kwargs: dict) -> None:
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
    def data():
        doc = "Access to raw data frame object."
        def fget(self):
            """Get dataframe"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set dataframe"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete dataframe"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    data = property(**data())

    ################## Column attributes ###################
    def ncol():
        doc = "Number of columns"
        def fget(self):
            """Get number of columns"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set number of columns"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete number of columns"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    ncol = property(**ncol())

    def col_axis():
        doc = "Column axis index"
        def fget(self):
            """Get column axis index"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set column axis index"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete column axis index"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    col_axis = property(**col_axis())

    def col_dtype():
        doc = "Column data types."
        def fget(self):
            """Get column data types"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set column data types"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete column data types"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    col_dtype = property(**col_dtype())

    def col_name():
        doc = "Column names."
        def fget(self):
            """Get column names"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set column names"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete column names"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    col_name = property(**col_name())

    def col_grp():
        doc = "Column groups used for classifying variables"
        def fget(self):
            """Get column groups"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set column groups"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete column groups"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    col_grp = property(**col_grp())

    #################### Row attributes ####################
    def nrow():
        doc = "Number of rows"
        def fget(self):
            """Get number of rows"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set number of rows"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete number of rows"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    nrow = property(**nrow())

    def row_axis():
        doc = "Row axis index"
        def fget(self):
            """Get row axis index"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set row axis index"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete row axis index"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    row_axis = property(**row_axis())

    def row_name():
        doc = "Row names."
        def fget(self):
            """Get row names"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set row names"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete row names"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    row_name = property(**row_name())

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
