import numpy
import pytest

from pybropt.popgen.gmat import DenseGenotypeMatrix
from pybropt.popgen.gmat import is_DenseGenotypeMatrix

@pytest.fixture
def mat_int8():
    a = numpy.int8([
        [0, 1, 1, 0, 0, 1, 0, 0],
        [1, 1, 0, 1, 1, 1, 1, 0],
        [1, 0, 0, 0, 0, 1, 1, 0],
        [1, 0, 0, 0, 1, 1, 1, 0],
        [0, 1, 1, 0, 1, 0, 1, 0],
        [0, 0, 1, 1, 0, 1, 0, 1],
        [0, 0, 0, 0, 0, 1, 0, 1],
        [1, 1, 1, 0, 1, 0, 1, 0]
    ])
    yield a

@pytest.fixture
def dgmat(mat_int8):
    yield DenseGenotypeMatrix(mat_int8)

def test_is_DenseGenotypeMatrix(dgmat):
    assert is_DenseGenotypeMatrix(dgmat)

# def test_math_operators(dgmat, mat_int8):
#     ops = [
#         "__add__","__sub__","__mul__",
#         "__truediv__","__floordiv__","__mod__","__pow__",
#         "__lshift__","__rshift__","__and__","__xor__","__or__",
#         "__radd__","__rsub__","__rmul__",
#         "__rtruediv__","__rfloordiv__","__rmod__",
#         "__rlshift__","__rrshift__","__rand__","__rxor__","__ror__",
#         "__iadd__","__isub__","__imul__","__imatmul__",
#         "__itruediv__","__ifloordiv__","__imod__","__ipow__",
#         "__ilshift__","__irshift__","__iand__","__ixor__","__ior__"
#     ]
#     ops_special = [
#         "__matmul__","__divmod__","__rmatmul__","__rdivmod__",
#     ]
#
#     for op in ops:
#         a = getattr(dgmat, op)(2)
#         b = getattr(mat_int8, op)(2)
#         assert numpy.all(a == b)
