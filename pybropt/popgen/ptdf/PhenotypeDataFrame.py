from pybropt.core.df import DataFrame

class PhenotypeDataFrame(DataFrame):
    """Abstract class for phenotype dataframe objects."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(PhenotypeDataFrame, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def analysis_type():
        doc = "Analysis variable type array."
        def fget(self):
            """Get analysis variable type array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set analysis variable type array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete analysis variable type array"""
            raise NotImplementedError("method is abstract")
        return locals()
    analysis_type = property(**analysis_type())

    def analysis_effect():
        doc = "Analysis variable effect type {'response','fixed','random',None} array."
        def fget(self):
            """Get analysis variable effect type array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set analysis variable effect type array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete analysis variable effect type array"""
            raise NotImplementedError("method is abstract")
        return locals()
    analysis_effect = property(**analysis_effect())



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
