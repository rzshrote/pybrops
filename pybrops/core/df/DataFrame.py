"""
Module defining DataFrame interfaces and associated error checking routines.
"""

class DataFrame:
    """
    An abstract class for data frame wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Data frame row and column access
        2) Data extraction from a data frame
        3) Data frame conversion
    """

    ########################## Special Object Methods ##########################
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

    ############################ Object Properties #############################
    @property
    def data(self) -> object:
        """Access to raw data frame object."""
        raise NotImplementedError("property is abstract")
    @data.setter
    def data(self, value: object) -> None:
        """Set dataframe"""
        raise NotImplementedError("property is abstract")

    ################## Column attributes ###################
    @property
    def ncol(self) -> object:
        """Number of columns."""
        raise NotImplementedError("property is abstract")
    @ncol.setter
    def ncol(self, value: object) -> None:
        """Set number of columns"""
        raise NotImplementedError("property is abstract")

    @property
    def col_axis(self) -> int:
        """Column axis index."""
        raise NotImplementedError("property is abstract")
    @col_axis.setter
    def col_axis(self, value: int) -> None:
        """Set column axis index"""
        raise NotImplementedError("property is abstract")

    @property
    def col_dtype(self) -> object:
        """Column data types."""
        raise NotImplementedError("property is abstract")
    @col_dtype.setter
    def col_dtype(self, value: object) -> None:
        """Set column data types"""
        raise NotImplementedError("property is abstract")
    
    @property
    def col_name(self) -> object:
        """Column names."""
        raise NotImplementedError("property is abstract")
    @col_name.setter
    def col_name(self, value: object) -> None:
        """Set column names"""
        raise NotImplementedError("property is abstract")
    
    @property
    def col_grp(self) -> object:
        """Column groups used for classifying variables."""
        raise NotImplementedError("property is abstract")
    @col_grp.setter
    def col_grp(self, value: object) -> None:
        """Set column groups"""
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

    @property
    def row_axis(self) -> int:
        """Row axis index."""
        raise NotImplementedError("property is abstract")
    @row_axis.setter
    def row_axis(self, value: int) -> None:
        """Set row axis index"""
        raise NotImplementedError("property is abstract")

    @property
    def row_name(self) -> object:
        """Row names."""
        raise NotImplementedError("property is abstract")
    @row_name.setter
    def row_name(self, value: object) -> None:
        """Set row names"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################
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



################################## Utilities ###################################
def check_is_DataFrame(v: object, vname: str) -> None:
    """
    Check if object is of type DataFrame. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DataFrame):
        raise TypeError("variable '{0}' must be a DataFrame".format(vname))
