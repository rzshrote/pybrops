"""
Module providing a dense coancestry matrix implementation using the VanRaden
method and associated error checking routines.
"""

from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix

# TODO: implement VanRaden methods
class DenseVanRadenCoancestryMatrix(DenseCoancestryMatrix):
    """
    A concrete class for a dense coancestry matrix calculated using the VanRaden
    method. Coancestry matrices are square.

    The purpose of this concrete class is to implement functionality for:
        1) Dense coancestry matrix value calculation.
        2) Dense coancestry matrix value access.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the concrete class DenseVanRadenCoancestryMatrix.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseVanRadenCoancestryMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
