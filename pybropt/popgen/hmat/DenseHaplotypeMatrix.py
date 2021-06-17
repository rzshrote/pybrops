from . import HaplotypeMatrix
from pybropt.core.error import check_is_ndarray
from pybropt.core.mat import DenseMutableMatrix

class DenseHaplotypeMatrix(DenseMutableMatrix,HaplotypeMatrix):
    """docstring for DenseHaplotypeMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, **kwargs):
        super(DenseHaplotypeMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Haplotype Data Properites ##############
    def mat():
        doc = "The mat property."
        def fget(self):
            return self._mat
        def fset(self, value):
            # The only assumption is that mat is a numpy.ndarray matrix.
            # Let the user decide whether to overwrite error checks.
            check_is_ndarray(value, "mat")
            self._mat = value
        def fdel(self):
            del self._mat
        return locals()
    mat = property(**mat())



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseHaplotypeMatrix(v):
    return isinstance(v, DenseHaplotypeMatrix)

def check_is_DenseHaplotypeMatrix(v, varname):
    if not isinstance(v, DenseHaplotypeMatrix):
        raise TypeError("'%s' must be a DenseHaplotypeMatrix." % varname)

def cond_check_is_DenseHaplotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DenseHaplotypeMatrix(v, varname)
