from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix

class AdditiveGenicVarianceMatrix(GenicVarianceMatrix):
    """docstring for AdditiveGenicVarianceMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class AdditiveGenicVarianceMatrix.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments. Used for cooperative inheritance.
            Dictionary passing unused arguments to the parent class constructor.
        """
        super(AdditiveGenicVarianceMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    @classmethod
    def from_algmod(cls, algmod, pgmat, nprogeny, mem):
        """
        Estimate genetic variances from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genetic
            variance.
        mem : int
            Memory chunk size to use during matrix operations.

        Returns
        -------
        out : GenicVarianceMatrix
            A matrix of additive genic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_AdditiveGenicVarianceMatrix(obj):
    """
    Determine whether an object is a ``AdditiveGenicVarianceMatrix``.

    Parameters
    ----------
    obj : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``AdditiveGenicVarianceMatrix`` object instance.
    """
    return isinstance(obj, AdditiveGenicVarianceMatrix)

def check_is_AdditiveGenicVarianceMatrix(obj, objname):
    """
    Check if object is of type ``AdditiveGenicVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(obj, AdditiveGenicVarianceMatrix):
        raise TypeError("'{0}' must be a AdditiveGenicVarianceMatrix".format(objname))

def cond_check_is_AdditiveGenicVarianceMatrix(obj, objname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type ``AdditiveGenicVarianceMatrix``. Otherwise raise
    ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    cond : function
        A function returning ``True`` or ``False`` for whether to test if ``obj``
        is a ``AdditiveGenicVarianceMatrix``.
    """
    if cond(obj) and not isinstance(obj, AdditiveGenicVarianceMatrix):
        raise TypeError("'{0}' must be a AdditiveGenicVarianceMatrix".format(objname))