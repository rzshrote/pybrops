from pybropt.core.mat.SquareTaxaMatrix import SquareTaxaMatrix

class GeneticVarianceMatrix(SquareTaxaMatrix):
    """docstring for GeneticVarianceMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class GeneticVarianceMatrix.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments. Used for cooperative inheritance.
            Dictionary passing unused arguments to the parent class constructor.
        """
        super(GeneticVarianceMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    @classmethod
    def from_gmod(cls, gmod, pgmat, ncross, nprogeny, s):
        """
        Estimate genetic variances from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : int
            Number of cross patterns to simulate for genetic variance
            estimation.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genetic
            variance.
        s : int
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.

        Returns
        -------
        out : GeneticVarianceMatrix
            A matrix of genetic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GeneticVarianceMatrix(obj):
    """
    Determine whether an object is a ``GeneticVarianceMatrix``.

    Parameters
    ----------
    obj : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``GeneticVarianceMatrix`` object instance.
    """
    return isinstance(obj, GeneticVarianceMatrix)

def check_is_GeneticVarianceMatrix(obj, objname):
    """
    Check if object is of type ``GeneticVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(obj, GeneticVarianceMatrix):
        raise TypeError("'{0}' must be a GeneticVarianceMatrix".format(objname))

def cond_check_is_GeneticVarianceMatrix(obj, objname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type ``GeneticVarianceMatrix``. Otherwise raise
    ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    cond : function
        A function returning ``True`` or ``False`` for whether to test if ``obj``
        is a ``GeneticVarianceMatrix``.
    """
    if cond(obj) and not isinstance(obj, GeneticVarianceMatrix):
        raise TypeError("'{0}' must be a GeneticVarianceMatrix".format(objname))
