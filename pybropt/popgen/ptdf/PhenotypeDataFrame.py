from pybropt.core.df.DataFrame import DataFrame

class PhenotypeDataFrame(DataFrame):
    """Abstract class for phenotype dataframe objects."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class PhenotypeDataFrame.
        """
        super(PhenotypeDataFrame, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    # TODO: maybe eliminate these. it seems like this information should go in other modules
    def col_analysis_type():
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
    col_analysis_type = property(**col_analysis_type())

    # TODO: maybe eliminate these. it seems like this information should go in other modules
    def col_analysis_effect():
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
    col_analysis_effect = property(**col_analysis_effect())



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhenotypeDataFrame(v):
    """
    Determine whether an object is a PhenotypeDataFrame.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a PhenotypeDataFrame object instance.
    """
    return isinstance(v, PhenotypeDataFrame)

def check_is_PhenotypeDataFrame(v, vname):
    """
    Check if object is of type PhenotypeDataFrame. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, PhenotypeDataFrame):
        raise TypeError("variable '{0}' must be a PhenotypeDataFrame".format(vname))

def cond_check_is_PhenotypeDataFrame(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type PhenotypeDataFrame. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a PhenotypeDataFrame.
    """
    if cond(v):
        check_is_PhenotypeDataFrame(v, vname)
