from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix

# TODO: implement VanRaden methods
class DenseVanRadenCoancestryMatrix(DenseCoancestryMatrix):
    """docstring for DenseVanRadenCoancestryMatrix."""

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
