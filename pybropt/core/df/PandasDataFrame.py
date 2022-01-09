from . import DataFrame

from pybropt.core.error import check_is_pandas_df

class PandasDataFrame(DataFrame):
    """docstring for PandasDataFrame."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, df, col_name = None, col_ctype = None, col_dtype = None, **kwargs):
        super(PandasDataFrame, self).__init__(**kwargs)
        self.df = df
        if col_name is not None:
            self.col_name = col_name
        self.col_ctype = col_ctype
        if col_dtype is not None:
            self.col_dtype = col_dtype

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def df():
        doc = "Access to raw dataframe."
        def fget(self):
            """Get dataframe"""
            return self._df
        def fset(self, value):
            """Set dataframe"""
            check_is_pandas_df(value, "df")
            self._df = value
        def fdel(self):
            """Delete dataframe"""
            del self._df
        return locals()
    df = property(**df())

    ################## Column attributes ###################
    def ncol():
        doc = "Number of columns"
        def fget(self):
            """Get number of columns"""
            return len(self._df.columns)
        def fset(self, value):
            """Set number of columns"""
            error_readonly("ncol")
        def fdel(self):
            """Delete number of columns"""
            error_readonly("ncol")
        return locals()
    ncol = property(**ncol())

    def col_axis():
        doc = "Column axis index"
        def fget(self):
            """Get column axis index"""
            return 1
        def fset(self, value):
            """Set column axis index"""
            error_readonly("col_axis")
        def fdel(self):
            """Delete column axis index"""
            error_readonly("col_axis")
        return locals()
    col_axis = property(**col_axis())

    def col_dtype():
        doc = "Column data types."
        def fget(self):
            """Get column data types"""
            return self._df.dtypes.values
        def fset(self, value):
            """Set column data types"""
            if isinstance(value, (list,tuple,numpy.ndarray)):
                check_len(value, "col_dtype", self.ncol)     # check input length
                names = self.col_name                        # get column names
                dtypes = self.col_dtype                      # get column dtypes
                # construct dict of different types
                value = dict((a,v) for a,b,v in zip(names,dtypes,value) if b != v)
            if isinstance(value, dict):
                self._df = self._df.astype(value)           # convert data types
            else:
                raise TypeError("unknown type: available types are numpy.ndarray, list, tuple, and dict")
        def fdel(self):
            """Delete column data types"""
            raise NotImplementedError("method is abstract")
        return locals()
    col_dtype = property(**col_dtype())

    def col_name():
        doc = "Column names."
        def fget(self):
            """Get column names"""
            return self._df.columns.values
        def fset(self, value):
            """Set column names"""
            self._df.columns = value
        def fdel(self):
            """Delete column names"""
            del self._df.columns
        return locals()
    col_name = property(**col_name())

    def col_ctype():
        doc = "Column types used for classifying variables"
        def fget(self):
            """Get column types"""
            return self._col_ctype
        def fset(self, value):
            """Set column types"""
            check_is_ndarray(value, "col_ctype")
            check_len(value, "col_ctype", self.ncol)
            self._col_ctype = value
        def fdel(self):
            """Delete column types"""
            del self._col_ctype
        return locals()
    col_ctype = property(**col_ctype())

    #################### Row attributes ####################
    def nrow():
        doc = "Number of rows"
        def fget(self):
            """Get number of rows"""
            return len(self._df.index)
        def fset(self, value):
            """Set number of rows"""
            error_readonly("nrow")
        def fdel(self):
            """Delete number of rows"""
            error_readonly("nrow")
        return locals()
    nrow = property(**nrow())

    def row_axis():
        doc = "Row axis index"
        def fget(self):
            """Get row axis index"""
            return 0
        def fset(self, value):
            """Set row axis index"""
            error_readonly("row_axis")
        def fdel(self):
            """Delete row axis index"""
            error_readonly("row_axis")
        return locals()
    row_axis = property(**row_axis())

    def row_name():
        doc = "Row names."
        def fget(self):
            """Get row names"""
            return self._df.index.values
        def fset(self, value):
            """Set row names"""
            self._df.index = value
        def fdel(self):
            """Delete row names"""
            del self._df.index
        return locals()
    row_name = property(**row_name())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def col_data(self, index = None, name = None, ctype = None, dtype = None, return_index = False, return_name = False, return_ctype = False, return_dtype = False, **kwargs):
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
        out_arr = [self.df.iloc[:,ix].values for ix in cix]

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

    def to_pandas_df(self, **kwargs):
        """
        Get dataframe as a pandas.DataFrame.
        """
        return self._df

    def to_dict(self, **kwargs):
        """
        Get dataframe as a dictionary of numpy.ndarray's.
        """
        df = self._df                                   # get pointer to pandas.DataFrame
        d = dict((e, df[e].values) for e in df.columns) # construct dictionary: col_name,array
        return d



################################################################################
################################## Utilities ###################################
################################################################################
def is_PandasDataFrame(v):
    return isinstance(v, PandasDataFrame)

def check_is_PandasDataFrame(v, vname):
    if not isinstance(v, PandasDataFrame):
        raise TypeError("variable '{0}' must be a PandasDataFrame".format(vname))

def cond_check_is_PandasDataFrame(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_PandasDataFrame(v, vname)
