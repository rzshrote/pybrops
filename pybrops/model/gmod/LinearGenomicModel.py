from pybrops.model.gmod.GenomicModel import GenomicModel

# TODO: change how linear models are defined so that they are inclusive of
#       dominance and epistatic genetic control architectures
class LinearGenomicModel(GenomicModel):
    """
    The LinearGenomicModel class represents a Multivariate Multiple Linear
    Regression model.

    A Multivariate Multiple Linear Regression model is defined as:

    .. math::
        Y = X \\Beta + ZU + e

    Where:

    - :math:`Y` is a matrix of response variables of shape ``(n,t)``.
    - :math:`X` is a matrix of fixed effect predictors of shape ``(n,q)``.
    - :math:`\\Beta` is a matrix of fixed effect regression coefficients of
      shape ``(q,t)``.
    - :math:`Z` is a matrix of random effect predictors of shape ``(n,p)``.
    - :math:`U` is a matrix of random effect regression coefficients of shape
      ``(p,t)``.
    - :math:`e` is a matrix of error terms of shape ``(n,t)``.

    Shape definitions:

    - ``n`` is the number of individuals
    - ``q`` is the number of fixed effect predictors (e.g. environments)
    - ``p`` is the number of random effect predictors (e.g. genomic markers)
    - ``t`` is the number of traits
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for LinearGenomicModel class.

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(LinearGenomicModel, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Linear Genomic Model Data ###############
    def beta():
        doc = "Fixed effect regression coefficients"
        def fget(self):
            """Get fixed effect regression coefficients"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set fixed effect regression coefficients"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete fixed effect regression coefficients"""
            raise NotImplementedError("method is abstract")
        return locals()
    beta = property(**beta())

    def u():
        doc = "Random effect regression coefficients"
        def fget(self):
            """Get random effect regression coefficients"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set random effect regression coefficients"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete random effect regression coefficients"""
            raise NotImplementedError("method is abstract")
        return locals()
    u = property(**u())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_LinearGenomicModel(v):
    """
    Determine whether an object is a LinearGenomicModel.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a LinearGenomicModel object instance.
    """
    return isinstance(v, LinearGenomicModel)

def check_is_LinearGenomicModel(v, vname):
    """
    Check if object is of type LinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, LinearGenomicModel):
        raise TypeError("variable '{0}' must be a LinearGenomicModel".format(vname))

def cond_check_is_LinearGenomicModel(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type LinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a LinearGenomicModel.
    """
    if cond(v):
        check_is_LinearGenomicModel(v, vname)