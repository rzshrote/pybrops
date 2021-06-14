from . import HaplotypeMatrix
from pybropt.core.error import check_is_ndarray

# TODO: inherit from DenseMutableMatrix
class DenseHaplotypeMatrix(HaplotypeMatrix):
    """docstring for DenseHaplotypeMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, **kwargs):
        super(DenseHaplotypeMatrix, self).__init__(
            mat = mat,
            **kwargs
        )
        self.mat = mat

    ############## Forward numeric operators ###############
    def __add__(self, value):
        return self._mat + value

    def __sub__(self, value):
        return self._mat - value

    def __mul__(self, value):
        return self._mat * value

    def __matmul__(self, value):
        return self._mat @ value

    def __truediv__(self, value):
        return self._mat / value

    def __floordiv__(self, value):
        return self._mat // value

    def __mod__(self, value):
        return self._mat % value

    def __divmod__(self, value):
        return divmod(self._mat, value)

    def __pow__(self, value):
        return self._mat ** value

    def __lshift__(self, value):
        return self._mat << value

    def __rshift__(self, value):
        return self._mat >> value

    def __and__(self, value):
        return self._mat & value

    def __xor__(self, value):
        return self._mat ^ value

    def __or__(self, value):
        return self._mat | value

    ############# Backwards numeric operators ##############
    def __radd__(self, value):
        return value + self._mat

    def __rsub__(self, value):
        return value - self._mat

    def __rmul__(self, value):
        return value * self._mat

    def __rmatmul__(self, value):
        return value @ self._mat

    def __rtruediv__(self, value):
        return value / self._mat

    def __rfloordiv__(self, value):
        return value // self._mat

    def __rmod__(self, value):
        return value % self._mat

    def __rdivmod__(self, value):
        return divmod(value, self._mat)

    def __rlshift__(self, value):
        return value << self._mat

    def __rrshift__(self, value):
        return value >> self._mat

    def __rand__(self, value):
        return value & self._mat

    def __rxor__(self, value):
        return value ^ self._mat

    def __ror__(self, value):
        return value | self._mat

    ############# Augmented numeric operators ##############
    def __iadd__(self, value):
        self._mat += value

    def __isub__(self, value):
        self._mat -= value

    def __imul__(self, value):
        self._mat *= value

    def __imatmul__(self, value):
        self._mat @= value

    def __itruediv__(self, value):
        self._mat /= value

    def __ifloordiv__(self, value):
        self._mat //= value

    def __imod__(self, value):
        self._mat %= value

    def __ipow__(self, value):
        self._mat **= value

    def __ilshift__(self, value):
        self._mat <<= value

    def __irshift__(self, value):
        self._mat >>= value

    def __iand__(self, value):
        self._mat &= value

    def __ixor__(self, value):
        self._mat ^= value

    def __ior__(self, value):
        self._mat |= value

    ################## Logical operators ###################
    def __lt__(self, value):
        return self._mat < value

    def __le__(self, value):
        return self._mat <= value

    def __eq__(self, value):
        return self._mat == value

    def __ne__(self, value):
        return self._mat != value

    def __gt__(self, value):
        return self._mat > value

    def __ge__(self, value):
        return self._mat >= value

    ################# Container operators ##################
    def __len__(self):
        return len(self._mat)

    def __getitem__(self, key):
        return self._mat[key]

    def __setitem__(self, key, value):
        self._mat[key] = value

    def __delitem__(self, key):
        del self._mat[key]

    def __iter__(self):
        return iter(self._mat)

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
