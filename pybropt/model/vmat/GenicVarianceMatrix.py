from pybropt.core.mat.SquareTaxaMatrix import SquareTaxaMatrix

class GenicVarianceMatrix(SquareTaxaMatrix):
    """docstring for GenicVarianceMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class GenicVarianceMatrix.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments. Used for cooperative inheritance.
            Dictionary passing unused arguments to the parent class constructor.
        """
        super(GenicVarianceMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    @classmethod
    def from_gmod(cls, gmod, pgmat, nprogeny):
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

        Returns
        -------
        out : GeneticVarianceMatrix
            A matrix of genic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenicVarianceMatrix(obj):
    """
    Determine whether an object is a ``GenicVarianceMatrix``.

    Parameters
    ----------
    obj : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``GenicVarianceMatrix`` object instance.
    """
    return isinstance(obj, GenicVarianceMatrix)

def check_is_GenicVarianceMatrix(obj, objname):
    """
    Check if object is of type ``GenicVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(obj, GenicVarianceMatrix):
        raise TypeError("'{0}' must be a GenicVarianceMatrix".format(objname))

def cond_check_is_GenicVarianceMatrix(obj, objname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type ``GenicVarianceMatrix``. Otherwise raise
    ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    cond : function
        A function returning ``True`` or ``False`` for whether to test if ``obj``
        is a ``GenicVarianceMatrix``.
    """
    if cond(obj) and not isinstance(obj, GenicVarianceMatrix):
        raise TypeError("'{0}' must be a GenicVarianceMatrix".format(objname))
