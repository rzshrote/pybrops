"""
Module defining interfaces and associated error checking routines for
optimization algorithms.
"""

class OptimizationAlgorithm:
    """
    An abstract class for optimization algorithms.

    The purpose of this abstract class is to provide functionality for:
        1) Optimization of objective functions.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class OptimizationAlgorithm.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(OptimizationAlgorithm, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def optimize(self, objfn, k, sspace, objfn_wt, **kwargs):
        """
        Optimize an objective function.

        Parameters
        ----------
        objfn : callable
            Objective function which to optimize.
        k : int
            Number of decision variables in the search space.
            A vector is formed as sspace^k
        sspace : numpy.ndarray
            Search space that the OptimizationAlgorithm searches in.
        objfn_wt : numpy.ndarray
            Weight(s) applied to output(s) from the objfn.

        Returns
        -------
        out : tuple
            A tuple of length 3 (soln, decn, misc)
        """
        raise NotImplementedError("method is abstract")

    # def optimize_vec(fn):
    #     raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_OptimizationAlgorithm(v):
    """
    Determine whether an object is a OptimizationAlgorithm.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a OptimizationAlgorithm object instance.
    """
    return isinstance(v, OptimizationAlgorithm)

def check_is_OptimizationAlgorithm(v, vname):
    """
    Check if object is of type OptimizationAlgorithm. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, OptimizationAlgorithm):
        raise TypeError("variable '{0}' must be a OptimizationAlgorithm".format(vname))

def cond_check_is_OptimizationAlgorithm(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type OptimizationAlgorithm. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a OptimizationAlgorithm.
    """
    if cond(v):
        check_is_OptimizationAlgorithm(v, vname)
