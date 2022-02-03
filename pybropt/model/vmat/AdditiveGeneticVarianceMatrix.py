from pybropt.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix

class AdditiveGeneticVarianceMatrix(GeneticVarianceMatrix):
    """docstring for AdditiveGeneticVarianceMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class AdditiveGeneticVarianceMatrix.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments. Used for cooperative inheritance.
            Dictionary passing unused arguments to the parent class constructor.
        """
        super(AdditiveGeneticVarianceMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    @classmethod
    def from_algmod(cls, algmod, pgmat, ncross, nprogeny, s, gmapfn, mem):
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
        gmapfn : GeneticMapFunction
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        mem : int
            Memory chunk size to use during matrix operations.

        Returns
        -------
        out : GeneticVarianceMatrix
            A matrix of additive genetic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_AdditiveGeneticVarianceMatrix(obj):
    """
    Determine whether an object is a ``AdditiveGeneticVarianceMatrix``.

    Parameters
    ----------
    obj : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``AdditiveGeneticVarianceMatrix`` object instance.
    """
    return isinstance(obj, AdditiveGeneticVarianceMatrix)

def check_is_AdditiveGeneticVarianceMatrix(obj, objname):
    """
    Check if object is of type ``AdditiveGeneticVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(obj, AdditiveGeneticVarianceMatrix):
        raise TypeError("'{0}' must be a AdditiveGeneticVarianceMatrix".format(objname))

def cond_check_is_AdditiveGeneticVarianceMatrix(obj, objname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type ``AdditiveGeneticVarianceMatrix``. Otherwise raise
    ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    cond : function
        A function returning ``True`` or ``False`` for whether to test if ``obj``
        is a ``AdditiveGeneticVarianceMatrix``.
    """
    if cond(obj) and not isinstance(obj, AdditiveGeneticVarianceMatrix):
        raise TypeError("'{0}' must be a AdditiveGeneticVarianceMatrix".format(objname))
