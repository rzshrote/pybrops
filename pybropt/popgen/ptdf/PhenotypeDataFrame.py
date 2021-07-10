class PhenotypeDataFrame:
    """Abstract class for phenotype dataframe objects."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(PhenotypeDataFrame, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def df():
        doc = "Access to raw data frame."
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

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def to_pandas_df(self, **kwargs):
        """
        Get dataframe as a pandas.DataFrame.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhenotypeDataFrame(v):
    return isinstance(v, PhenotypeDataFrame)

def check_is_PhenotypeDataFrame(v, vname):
    if not isinstance(v, PhenotypeDataFrame):
        raise TypeError("variable '{0}' must be a PhenotypeDataFrame".format(vname))

def cond_check_is_PhenotypeDataFrame(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_PhenotypeDataFrame(v, vname)
