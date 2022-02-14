from pybrops.model.vmat.DenseGeneticVarianceMatrix import DenseGeneticVarianceMatrix
from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix

# TODO: implement me
class DenseAdditiveGeneticVarianceMatrix(DenseGeneticVarianceMatrix,AdditiveGeneticVarianceMatrix):
    """docstring for DenseAdditiveGeneticVarianceMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, **kwargs):
        """
        Constructor for the concrete class DenseAdditiveGeneticVarianceMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            Array used to construct the object.
        taxa : numpy.ndarray
            Taxa names.
        taxa_grp : numpy.ndarray
            Taxa groupings.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseAdditiveGeneticVarianceMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )
    
    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    # from_gmod
    # this method should remain abstract; it depends on the cross structure
    # maybe in the future if generic cross structures are implemented,
    # then this method can become concrete

    # from_algmod
    # this method should remain abstract; it depends on the cross structure
    # maybe in the future if generic cross structures are implemented,
    # then this method can become concrete



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseAdditiveGeneticVarianceMatrix(obj):
    """
    Determine whether an object is a ``DenseAdditiveGeneticVarianceMatrix``.

    Parameters
    ----------
    obj : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``DenseAdditiveGeneticVarianceMatrix`` object instance.
    """
    return isinstance(obj, DenseAdditiveGeneticVarianceMatrix)

def check_is_DenseAdditiveGeneticVarianceMatrix(obj, objname):
    """
    Check if object is of type ``DenseAdditiveGeneticVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(obj, DenseAdditiveGeneticVarianceMatrix):
        raise TypeError("'{0}' must be a DenseAdditiveGeneticVarianceMatrix".format(objname))

def cond_check_is_DenseAdditiveGeneticVarianceMatrix(obj, objname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type ``DenseAdditiveGeneticVarianceMatrix``. Otherwise raise
    ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    cond : function
        A function returning ``True`` or ``False`` for whether to test if ``obj``
        is a ``DenseAdditiveGeneticVarianceMatrix``.
    """
    if cond(obj) and not isinstance(obj, DenseAdditiveGeneticVarianceMatrix):
        raise TypeError("'{0}' must be a DenseAdditiveGeneticVarianceMatrix".format(objname))
