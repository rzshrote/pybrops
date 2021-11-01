from pybropt.core.mat import TaxaMatrix

class CoancestryMatrix(TaxaMatrix):
    """docstring for CoancestryMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for CoancestryMatrix class.

        Parameters
        ----------
        **kwargs : **dict
            Additional keyword arguments for dependency injection.
        """
        super(CoancestryMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############### Genotype Data Properites ###############
    # gmat should be implemented in a sparse version of this matrix.

    ############## Coancestry Data Properites ##############
    # access using mat (inherited from Matrix)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################## Coancestry Methods ##################
    def coancestry(self, *args, **kwargs):
        """
        Retrieve the coancestry between individuals.

        Parameters
        ----------
        *args : *tuple
            A tuple of matrix indices to access the coancestry.
        **kwargs : **dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    @classmethod
    def from_gmat(cls, gmat, **kwargs):
        """
        Create a CoancestryMatrix from a GenotypeMatrix.

        Parameters
        ----------
        gmat : GenotypeMatrix
            Input genotype matrix from which to calculate coancestry.
        **kwargs : **dict
            Additional keyword arguments.
        """
        raise NotImplementedError("class method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_CoancestryMatrix(v):
    """
    Determine whether an object is a CoancestryMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a CoancestryMatrix object instance.
    """
    return isinstance(v, CoancestryMatrix)

def check_is_CoancestryMatrix(v, vname):
    """
    Check if object is of type CoancestryMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, CoancestryMatrix):
        raise TypeError("variable '{0}' must be a CoancestryMatrix".format(vname))

def cond_check_is_CoancestryMatrix(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type CoancestryMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        CoancestryMatrix.
    """
    if cond(v):
        check_is_CoancestryMatrix(v, vname)
