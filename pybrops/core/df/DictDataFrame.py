"""
Module implementing a Dictionary DataFrame and associated error checking routines.
"""

import copy
import numbers
from typing import Any, Union
import numpy
import pandas

from pybrops.core.df.DataFrame import DataFrame
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.core.error.error_value_python import check_keys_in_dict
from pybrops.core.error.error_value_python import check_keys_in_dict_all_type
from pybrops.core.error.error_value_python import check_len
from pybrops.core.error.error_value_python import check_values_in_dict_all_type
from pybrops.core.error.error_value_python import check_values_in_dict_equal_len
from pybrops.core.error.error_value_python import check_values_in_dict_len
from pybrops.core.error.error_attr_python import error_readonly

class DictDataFrame(DataFrame):
    """
    A concrete class for data frame objects utilizing Python dictionaries as a
    storage container.
    """

    ########################## Special Object Methods ##########################
    def __init__(self, data, col_grp = None, row_name = None, **kwargs: dict):
        """
        Constructor for the concrete class DictDataFrame.

        Parameters
        ----------
        data : dict
            Dictionary containing DataFrame's contents.
            Dictionary keys must be all strings.
            Dictionary values must be all numpy.ndarray's.
        col_grp : dict
            Column groups
        kwargs : dict
            Additional keyword arguments.
        """
        super(DictDataFrame, self).__init__()

        ### error checks and assignments (order dependent)
        self.data = data
        self.col_grp = col_grp
        self.row_name = row_name

    def __copy__(self):
        """
        Make a shallow copy of the the dataframe.

        Returns
        -------
        out : DataFrame
        """
        return self.__class__(
            data = copy.copy(self.data),
            col_grp = copy.copy(self.col_grp),
            row_name = copy.copy(self.row_name)
        )

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
        return self.__class__(
            data = copy.deepcopy(self.data),
            col_grp = copy.deepcopy(self.col_grp),
            row_name = copy.deepcopy(self.row_name)
        )

    ############################ Object Properties #############################
    @property
    def data(self) -> dict:
        """Access to raw data frame object."""
        return self._data
    @data.setter
    def data(self, value: dict) -> None:
        """Set dataframe"""
        check_is_dict(value, "data")
        check_keys_in_dict_all_type(value, "data", str)
        check_values_in_dict_all_type(value, "data", numpy.ndarray)
        check_values_in_dict_equal_len(value, "data")
        self._data = value

    ################## Column attributes ###################
    @property
    def ncol(self) -> int:
        """Number of columns."""
        return len(self._data)
    @ncol.setter
    def ncol(self, value: int) -> None:
        """Set number of columns"""
        error_readonly("ncol")

    @property
    def col_axis(self) -> int:
        """Column axis index."""
        return 1
    @col_axis.setter
    def col_axis(self, value: int) -> None:
        """Set column axis index"""
        error_readonly("col_axis")

    @property
    def col_dtype(self) -> numpy.ndarray:
        """Column data types."""
        return numpy.object_(list(e.dtype for e in self._data.values()))
    @col_dtype.setter
    def col_dtype(self, value: numpy.ndarray) -> None:
        """Set column data types"""
        if isinstance(value, (list,tuple,numpy.ndarray)):
            check_len(value, "col_dtype", self.ncol)    # check input length
            names = self.col_name                       # get column names
            dtypes = value                              # get column dtypes
            value = dict(zip(names,dtypes))             # construct dict of dtypes
        if isinstance(value, dict):
            names = self.col_name
            for k in value.keys():
                if isinstance(k, numbers.Integral):     # if key is an integer
                    if (k >= 0) and (k < len(names)):   # if integer range is acceptable
                        new_key = names[k]              # get new key
                        value[new_key] = value.pop(k)   # rename old key
                    else:                               # otherwise integer range is unacceptable
                        value.pop(k)                    # remove key from dictionary
                elif isinstance(k, str):                # if key is a string
                    if k not in names:                  # if key is not in column names
                        value.pop(k)                    # remove key from dictionary
                else:                                   # if key type not recognized
                    value.pop(k)                        # remove key from dictionary
            check_values_in_dict_all_type(value, "col_dtype", (str,numpy.dtype))
        else:
            raise TypeError("unsupported type: supported types are dict, list, tuple, numpy.ndarray")
        for k,v in value.items():                   # for each key, value pair
            if self._data[k].dtype != v:              # if dtypes do not match up
                self._data[k] = self._data[k].astype(v) # convert data types

    @property
    def col_name(self) -> numpy.ndarray:
        """Column names."""
        return numpy.object_(list(self._data.keys()))
    @col_name.setter
    def col_name(self, value: numpy.ndarray) -> None:
        """Set column names"""
        if isinstance(value, (list,tuple,numpy.ndarray)):
            check_len(value, "col_name", self.ncol)     # check input length
            names = self.col_name                       # get column names
            new_names = value                           # get new column names
            value = dict(zip(names,new_names))          # construct dict of names
        if isinstance(value, dict):
            names = self.col_name
            for k in value.keys():
                if isinstance(k, numbers.Integral):     # if key is an integer
                    if (k >= 0) and (k < len(names)):   # if integer range is acceptable
                        new_key = names[k]              # get new key
                        value[new_key] = value.pop(k)   # rename old key
                    else:                               # otherwise integer range is unacceptable
                        value.pop(k)                    # remove key from dictionary
                elif isinstance(k, str):                # if key is a string
                    if k not in names:                  # if key is not in column names
                        value.pop(k)                    # remove key from dictionary
                else:                                   # if key type not recognized
                    value.pop(k)                        # remove key from dictionary
            check_values_in_dict_all_type(value, "col_name", (str,numpy.dtype))
        else:
            raise TypeError("unsupported type: supported types are dict, list, tuple, numpy.ndarray")
        for k,v in value.items():                       # for each key, value pair
            self._data[v] = self._data.pop(k)           # rename old key

    @property
    def col_grp(self) -> Union[numpy.ndarray,None]:
        """Column groups used for classifying variables."""
        return self._col_grp
    @col_grp.setter
    def col_grp(self, value: Union[numpy.ndarray,None]) -> None:
        """Set column groups"""
        if isinstance(value, (list,tuple,numpy.ndarray)):   # if array_like
            check_len(value, "col_grp", self.ncol)          # check input length
            names = self.col_name                           # get column names
            grps = value                                    # get column grps
            value = dict(zip(names,grps))                   # construct dict of grps
        if value is None:                                   # if is None
            self._col_grp = value                           # set to None
        elif isinstance(value, dict):                       # if is dict
            check_keys_in_dict(value, "col_grp", *self._data.keys())
            check_values_in_dict_all_type(value, "col_grp", str)
            self._col_grp = {k: value[k] for k in self._data.keys()}
        else:
            raise TypeError("unsupported type: supported types are dict, list, tuple, numpy.ndarray")

    #################### Row attributes ####################
    @property
    def nrow(self) -> int:
        """Number of rows."""
        try:
            k0 = next(iter(self._data))
            return len(self._data[k0])
        except StopIteration:
            return 0
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
    def row_name(self) -> Union[numpy.ndarray,None]:
        """Row names."""
        if hasattr(self, "_row_name") and isinstance(self._row_name, dict):
            return self._row_name["row_name"]
        return None
    @row_name.setter
    def row_name(self, value: Union[numpy.ndarray,None]) -> None:
        """Set row names"""
        if isinstance(value, (list,tuple,numpy.ndarray)):
            check_len(value, "row_name", self.nrow)
            value = {"row_name": numpy.object_(value)}
        if value is None:
            self._row_name = value
        elif isinstance(value, dict):
            check_is_dict(value, "row_name")
            check_keys_in_dict(value, "row_name", "row_name")
            check_keys_in_dict_all_type(value, "row_name", str)
            check_values_in_dict_all_type(value, "row_name", numpy.ndarray)
            check_values_in_dict_len(value, "row_name", self.nrow)
            check_values_in_dict_equal_len(value, "row_name")
            self._row_name = value
        else:
            raise TypeError("unsupported type: supported types are dict, list, tuple, numpy.ndarray")

    ############################## Object Methods ##############################
    def col_data(
            self, 
            index = None, 
            name = None, 
            grp = None, 
            dtype = None, 
            return_index = False, 
            return_name = False, 
            return_grp = False, 
            return_dtype = False, 
            **kwargs: dict
        ):
        """
        Get a column's (or columns') data from the dataframe.

        Parameters
        ----------
        index : int, None
            Integer index of the column to get.
        name : str, None
            Name of the column to get.
        grp : str, None
            Classification type of the column to get.
        dtype : str, numpy.dtype, None
            Data type of the column to get.
        return_index : boolean, default = False
            Whether to return the column index along with the column data.
        return_name : boolean, default = False
            Whether to return the column name along with the column data.
        return_grp : boolean, default = False
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
        if grp is not None:
            mask = mask & (self.col_grp == grp)
        if dtype is not None:
            mask = mask & (self.col_dtype == dtype)
        if mask is True:
            mask = False

        # get index of column equal to the provided name
        cix = numpy.flatnonzero(mask)

        # process errors if key is not found or multiple keys are found
        if numpy.sum(mask) == 0:
            err =  "DictDataFrame does not contain a column with matches for: \n"
            err += "    name = {0}".format(name)
            err += "    grp = {0}".format(grp)
            err += "    dtype = {0}".format(dtype)
            raise KeyError(err)

        # extract column arrays as a numpy.ndarray
        out_arr = [e for e,m in zip(self.data.values(),mask) if m]

        # construct extra output list
        out_extra = []

        # add values to extra output
        if return_index:
            out_extra.append(cix)
        if return_name:
            out_extra.append(self.col_name[cix])
        if return_grp:
            out_extra.append(self.col_grp[cix])
        if return_dtype:
            out_extra.append(self.col_dtype[cix])

        # construct output object (numpy.ndarray or tuple)
        out = out_arr if len(out_extra) == 0 else (out_arr, *out_extra)

        return out

    def to_pandas_df(self, **kwargs: dict):
        """
        Get dataframe as a pandas.DataFrame.

        Returns
        -------
        out : pandas.DataFrame
            DataFrame as a pandas.DataFrame.
        """
        out = pandas.DataFrame(
            data = self.data,
            index = self.row_name
        )
        return out

    def to_dict(self, **kwargs: dict):
        """
        Get dataframe as a dictionary of numpy.ndarray's.

        Returns
        -------
        out : dict
            DataFrame as a dictionary of numpy.ndarray's.
        """
        return self._data



################################## Utilities ###################################
def check_is_DictDataFrame(v: object, vname: str) -> None:
    """
    Check if object is of type DictDataFrame. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DictDataFrame):
        raise TypeError("variable '{0}' must be a DictDataFrame".format(vname))
