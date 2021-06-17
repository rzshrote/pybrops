from . import CoancestryMatrix

from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_is_2d
from pybropt.core.error import check_ndarray_is_square
from pybropt.core.mat import DenseMutableMatrix

class DenseCoancestryMatrix(DenseMutableMatrix,CoancestryMatrix):
    """docstring for DenseCoancestryMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, **kwargs):
        super(DenseCoancestryMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Coancestry Data Properites ##############
    def mat():
        doc = "The mat property."
        def fget(self):
            return self._mat
        def fset(self, value):
            check_is_ndarray(value, "mat")          # must be numpy.ndarray
            check_ndarray_is_2d(value, "mat")       # must be 2D matrix
            check_ndarray_is_square(value, "mat")   # must be square
            self._mat = value
        def fdel(self):
            del self._mat
        return locals()
    mat = property(**mat())



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseCoancestryMatrix(v):
    return isinstance(v, DenseCoancestryMatrix)

def check_is_DenseCoancestryMatrix(v, varname):
    if not isinstance(v, DenseCoancestryMatrix):
        raise TypeError("'%s' must be a DenseCoancestryMatrix." % varname)

def cond_check_is_DenseCoancestryMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DenseCoancestryMatrix(v, varname)
