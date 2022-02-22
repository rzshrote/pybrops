"""
Module implementing matrix routines and associated error checking routines
for dense breeding value matrices estimated from phenotypic data.
"""

from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix

# TODO: add standard errors for this class; this could be used for two-stage estimation
class DenseEstimatedBreedingValueMatrix(DenseBreedingValueMatrix):
    """docstring for DenseEstimatedBreedingValueMatrix."""

    def __init__(self, mat, location, scale, taxa = None, taxa_grp = None, trait = None, **kwargs):
        """
        Constructor for the concrete class DenseEstimatedBreedingValueMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            A float64 matrix of breeding values of shape ``(n,t)``.
        location : numpy.ndarray
            A numpy.ndarray of shape ``(t,)`` containing breeding value locations.
        scale : numpy.ndarray
            A numpy.ndarray of shape ``(t,)`` containing breeding value scales.
        taxa : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.
        taxa_grp : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.
        trait : numpy.ndarray, None
            A numpy.ndarray of shape ``(t,)`` containing trait names.
            If ``None``, do not store any trait name information.
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(DenseEstimatedBreedingValueMatrix, self).__init__(
            mat = mat,
            location = location,
            scale = scale,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseEstimatedBreedingValueMatrix(v):
    """
    Determine whether an object is a DenseEstimatedBreedingValueMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseEstimatedBreedingValueMatrix object instance.
    """
    return isinstance(v, DenseEstimatedBreedingValueMatrix)

def check_is_DenseEstimatedBreedingValueMatrix(v, vname):
    """
    Check if object is of type DenseEstimatedBreedingValueMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseEstimatedBreedingValueMatrix):
        raise TypeError("variable '{0}' must be a DenseEstimatedBreedingValueMatrix".format(vname))

def cond_check_is_DenseEstimatedBreedingValueMatrix(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type DenseEstimatedBreedingValueMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        DenseEstimatedBreedingValueMatrix.
    """
    if cond(v):
        check_is_DenseEstimatedBreedingValueMatrix(v, vname)
