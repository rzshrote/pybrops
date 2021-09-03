class DataFrame:
    """docstring for DataFrame."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(DataFrame, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def df():
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
        return locals()
    df = property(**df())

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
        return locals()
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
        return locals()
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
        return locals()
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
        return locals()
    col_name = property(**col_name())

    def col_ctype():
        doc = "Column types used for classifying variables"
        def fget(self):
            """Get column types"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set column types"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete column types"""
            raise NotImplementedError("method is abstract")
        return locals()
    col_ctype = property(**col_ctype())

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
        return locals()
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
        return locals()
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
        return locals()
    row_name = property(**row_name())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def col_data(self, **kwargs):
        """
        Get a column's data from the dataframe.
        """
        raise NotImplementedError("method is abstract")

    def to_pandas_df(self, **kwargs):
        """
        Get dataframe as a pandas.DataFrame.
        """
        raise NotImplementedError("method is abstract")

    def to_dict(self, **kwargs):
        """
        Get dataframe as a dictionary of numpy.ndarray's.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_DataFrame(v):
    return isinstance(v, DataFrame)

def check_is_DataFrame(v, vname):
    if not isinstance(v, DataFrame):
        raise TypeError("variable '{0}' must be a DataFrame".format(vname))

def cond_check_is_DataFrame(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DataFrame(v, vname)
