from . import HaplotypeMatrix
from pybropt.core.mat import VariantMatrix
from pybropt.core.mat import TaxaMatrix
from pybropt.core.mat import SortableMatrix

class HaplotypeVariantMatrix(HaplotypeMatrix,VariantMatrix,TaxaMatrix,SortableMatrix):
    """docstring for HaplotypeVariantMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(HaplotypeVariantMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
