from pybropt.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix

# TODO: implement VanRaden methods
class DenseVanRadenCoancestryMatrix(DenseCoancestryMatrix):
    """docstring for DenseVanRadenCoancestryMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(DenseVanRadenCoancestryMatrix, self).__init__(
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
