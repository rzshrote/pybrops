"""
Module implementing a Pandas DataFrame and associated error checking routines.
"""

from typing import Any

import numpy
import pandas
from pybrops.core.df.DataFrame import DataFrame
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_python import check_len

class PandasDataFrame(DataFrame):
    """
    A concrete class for data frame objects utilizing Pandas DataFrames as a
    storage container.
    """

    ########################## Special Object Methods ##########################
    def __init__(self, data, col_name = None, col_ctype = None, col_dtype = None, **kwargs: dict):
        """
        Constructor for the concrete class PandasDataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
        """
        super(PandasDataFrame, self).__init__(**kwargs)
        self.data = data
        if col_name is not None:
            self.col_name = col_name
        self.col_ctype = col_ctype
        if col_dtype is not None:
            self.col_dtype = col_dtype

    ############################ Object Properties #############################
    @property
    def data(self) -> pandas.DataFrame:
        """Access to raw dataframe."""
        return self._data
    @data.setter
    def data(self, value: pandas.DataFrame) -> None:
        """Set dataframe"""
        check_is_pandas_DataFrame(value, "df")
        self._data = value

    ################## Column attributes ###################
    @property
    def ncol(self) -> int:
        """Number of columns."""
        return len(self._data.columns)
    @ncol.setter
    def ncol(self, value: int) -> None:
        """Set number of columns"""
        error_readonly("ncol")

    @property
    def col_axis(self) -> int:
        """Column axis index."""
        return self._col_axis
    @col_axis.setter
    def col_axis(self, value: int) -> None:
        """Set column axis index"""
        error_readonly("col_axis")
    
    @property
    def col_dtype(self) -> numpy.ndarray:
        """Column data types."""
        return self._col_dtype
    @col_dtype.setter
    def col_dtype(self, value: numpy.ndarray) -> None:
        """Set column data types"""
        if isinstance(value, (list,tuple,numpy.ndarray)):
            check_len(value, "col_dtype", self.ncol)     # check input length
            names = self.col_name                        # get column names
            dtypes = self.col_dtype                      # get column dtypes
            # construct dict of different types
            value = dict((a,v) for a,b,v in zip(names,dtypes,value) if b != v)
        if isinstance(value, dict):
            self._data = self._data.astype(value)           # convert data types
        else:
            raise TypeError("unknown type: available types are numpy.ndarray, list, tuple, and dict")
    
    @property
    def col_name(self) -> numpy.ndarray:
        """Column names."""
        return self._col_name
    @col_name.setter
    def col_name(self, value: numpy.ndarray) -> None:
        """Set column names"""
        self._data.columns = value

    @property
    def col_ctype(self) -> numpy.ndarray:
        """Column types used for classifying variables."""
        return self._col_ctype
    @col_ctype.setter
    def col_ctype(self, value: numpy.ndarray) -> None:
        """Set column types"""
        check_is_ndarray(value, "col_ctype")
        check_len(value, "col_ctype", self.ncol)
        self._col_ctype = value

    #################### Row attributes ####################
    @property
    def nrow(self) -> int:
        """Number of rows."""
        return len(self._data.index)
    @nrow.setter
    def nrow(self, value: int) -> None:
        """Set number of rows"""
        error_readonly("nrow")

    @property
    def row_axis(self) -> int:
        """Row axis index."""
        return 0
    @row_axis.setter
    def row_axis(self, value: int) -> None:
        """Set row axis index"""
        error_readonly("row_axis")

    @property
    def row_name(self) -> numpy.ndarray:
        """Row names."""
        return self._data.index.values
    @row_name.setter
    def row_name(self, value: numpy.ndarray) -> None:
        """Set row names"""
        self._data.index = value

    ############################## Object Methods ##############################
    def col_data(self, index = None, name = None, ctype = None, dtype = None, return_index = False, return_name = False, return_ctype = False, return_dtype = False, **kwargs: dict):
        """
        Get a column's (or columns') data from the dataframe.

        Parameters
        ----------
        index : int, None
            Integer index of the column to get.
        name : str, None
            Name of the column to get.
        ctype : str, None
            Classification type of the column to get.
        dtype : str, numpy.dtype, None
            Data type of the column to get.
        return_index : boolean, default = False
            Whether to return the column index along with the column data.
        return_name : boolean, default = False
            Whether to return the column name along with the column data.
        return_ctype : boolean, default = False
            Whether to return the column type along with the column data.
        return_dtype : boolean, default = False
            Whether to return the column dtype along with the column data.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray, tuple, None
            If no other returns are specified, return a numpy.ndarray for the
            column. Otherwise return a tuple of values corresponding to (in
            order of the function API) the requested data.
        """
        # construct a selection mask for column names
        mask = True
        if index is not None:
            tmp = numpy.empty(self.ncol, dtype = 'bool')
            tmp[index] = True
            mask = mask & tmp
        if name is not None:
            mask = mask & (self.col_name == name)
        if ctype is not None:
            mask = mask & (self.col_ctype == ctype)
        if dtype is not None:
            mask = mask & (self.col_dtype == dtype)
        if mask is True:
            mask = False

        # get index of column equal to the provided name
        cix = numpy.flatnonzero(self.col_name == name)

        # get number of indices found
        nix = len(cix)

        # process errors if key is not found or multiple keys are found
        if len(cix) == 0:
            err =  "PandasDataFrame does not contain a column with matches for: \n"
            err += "    name = {0}".format(name)
            err += "    ctype = {0}".format(ctype)
            err += "    dtype = {0}".format(dtype)
            raise KeyError(err)

        # extract column arrays as a numpy.ndarray
        out_arr = [self.data.iloc[:,ix].values for ix in cix]

        # construct extra output list
        out_extra = []

        # add values to extra output
        if return_index:
            out_extra.append(cix)
        if return_name:
            out_extra.append(self.col_name[cix])
        if return_ctype:
            out_extra.append(self.col_ctype[cix])
        if return_dtype:
            out_extra.append(self.col_dtype[cix])

        # construct output object (numpy.ndarray or tuple)
        out = out_arr if len(out_extra) == 0 else (out_arr, *out_extra)

        return out

    def to_pandas_df(self, **kwargs: dict):
        """
        Get dataframe as a pandas.DataFrame.
        """
        return self._data

    def to_dict(self, **kwargs: dict):
        """
        Get dataframe as a dictionary of numpy.ndarray's.
        """
        df = self._data                                   # get pointer to pandas.DataFrame
        d = dict((e, df[e].values) for e in df.columns) # construct dictionary: col_name,array
        return d



################################## Utilities ###################################
def check_is_PandasDataFrame(v: object, vname: str) -> None:
    if not isinstance(v, PandasDataFrame):
        raise TypeError("variable '{0}' must be a PandasDataFrame".format(vname))
