from . import PhenotypeDataFrame

from pybropt.core.error import check_is_pandas_df

class PandasPhenotypeDataFrame(PhenotypeDataFrame):
    """docstring for PandasPhenotypeDataFrame."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, df, **kwargs):
        super(PandasPhenotypeDataFrame, self).__init__(**kwargs)
        self.df = df

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

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def to_pandas_df(self, **kwargs):
        """
        Get dataframe as a pandas.DataFrame.
        """
        return self._df



################################################################################
################################## Utilities ###################################
################################################################################
def is_PandasPhenotypeDataFrame(v):
    return isinstance(v, PandasPhenotypeDataFrame)

def check_is_PandasPhenotypeDataFrame(v, vname):
    if not isinstance(v, PandasPhenotypeDataFrame):
        raise TypeError("variable '{0}' must be a PandasPhenotypeDataFrame".format(vname))

def cond_check_is_PandasPhenotypeDataFrame(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_PandasPhenotypeDataFrame(v, vname)
