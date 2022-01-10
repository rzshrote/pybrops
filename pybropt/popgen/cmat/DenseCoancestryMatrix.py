from pybropt.popgen.cmat.CoancestryMatrix import CoancestryMatrix

from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_is_2d
from pybropt.core.error import check_ndarray_is_square
from pybropt.core.mat import DenseTaxaMatrix

class DenseCoancestryMatrix(DenseTaxaMatrix,CoancestryMatrix):
    """docstring for DenseCoancestryMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, **kwargs):
        """
        Constructor for DenseCoancestryMatrix class.

        Parameters
        ----------
        mat : numpy.ndarray
            Coancestry matrix of shape (n,n)
            Where:
                n is the number of taxa.
        """
        # call constructor for DenseTaxaMatrix
        super(DenseCoancestryMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Coancestry Data Properites ##############
    def mat():
        doc = "The mat property."
        def fget(self):
            return self._mat
        def fset(self, value):
            check_is_ndarray(value, "mat")          # must be numpy.ndarray
            check_ndarray_is_2d(value, "mat")       # must be 2D matrix
            check_ndarray_is_square(value, "mat")   # must be square
            self._mat = value
        def fdel(self):
            del self._mat
        return locals()
    mat = property(**mat())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################## Coancestry Methods ##################
    def coancestry(self, *args, **kwargs):
        """
        Retrieve the coancestry between individuals.

        Parameters
        ----------
        args : tuple
            A tuple of matrix indices to access the coancestry.
        kwargs : dict
            Additional keyword arguments.
        """
        # index via numpy and return.
        return self._mat[args]

    # TODO: work on adjoin, delete, insert, select, concat, append, remove, etc.



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseCoancestryMatrix(v):
    """
    Determine whether an object is a DenseCoancestryMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseCoancestryMatrix object instance.
    """
    return isinstance(v, DenseCoancestryMatrix)

def check_is_DenseCoancestryMatrix(v, vname):
    """
    Check if object is of type DenseCoancestryMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseCoancestryMatrix):
        raise TypeError("variable '{0}' must be a DenseCoancestryMatrix".format(vname))

def cond_check_is_DenseCoancestryMatrix(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type DenseCoancestryMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        DenseCoancestryMatrix.
    """
    if cond(v):
        check_is_DenseCoancestryMatrix(v, vname)
