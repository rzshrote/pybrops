# our libraries
from . import GenomicModel

class LinearGenomicModel(GenomicModel):
    """
    The LinearGenomicModel class represents a Multivariate Multiple Linear
    Regression model.

    A Multivariate Multiple Linear Regression model is defined as:

        Y = 1u' + XB + e

        Where:
            Y is a matrix of response variables of shape (n, t)
            1 is a vector of ones of shape (n, 1)
            u is a vector of regression intercepts of shape (t, 1)
            X is a matrix of genotypes of shape (n, p)
            B is a matrix of regression coefficients of shape (p, t)
            e is a matrix of error terms of shape (n, t)

        Shape definitions:
            n = the number of individuals
            p = the number of markers
            t = the number of traits
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for LinearGenomicModel class.

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(LinearGenomicModel, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Linear Genomic Model Data ###############
    def mu():
        doc = "The mu property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    mu = property(**mu())

    def beta():
        doc = "The beta property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    beta = property(**beta())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
