import operator

from . import generic_test_numpy_operator

def generic_test_numpy_add(v1, v2, w1, w2):
    """
    Test:
        (v1 + v2) == (w1 + w2)
    """
    generic_test_numpy_operator(operator.add, (v1, v2), (w1, w2))

def generic_test_numpy_sub(v1, v2, w1, w2):
    """
    Test:
        (v1 - v2) == (w1 - w2)
    """
    generic_test_numpy_operator(operator.sub, (v1, v2), (w1, w2))

def generic_test_numpy_mul(v1, v2, w1, w2):
    """
    Test:
        (v1 * v2) == (w1 * w2)
    """
    generic_test_numpy_operator(operator.mul, (v1, v2), (w1, w2))

def generic_test_numpy_matmul(v1, v2, w1, w2):
    """
    Test:
        (v1 @ v2) == (w1 @ w2)
    """
    generic_test_numpy_operator(operator.matmul, (v1, v2), (w1, w2))

def generic_test_numpy_truediv(v1, v2, w1, w2):
    """
    Test:
        (v1 / v2) == (w1 / w2)
    """
    generic_test_numpy_operator(operator.truediv, (v1, v2), (w1, w2))

def generic_test_numpy_floordiv(v1, v2, w1, w2):
    """
    Test:
        (v1 // v2) == (w1 // w2)
    """
    generic_test_numpy_operator(operator.floordiv, (v1, v2), (w1, w2))

# TODO: divmod

def generic_test_numpy_pow(v1, v2, w1, w2):
    """
    Test:
        (v1 ** v2) == (w1 ** w2)
    """
    generic_test_numpy_operator(operator.pow, (v1, v2), (w1, w2))

def generic_test_numpy_lshift(v1, v2, w1, w2):
    """
    Test:
        (v1 << v2) == (w1 << w2)
    """
    generic_test_numpy_operator(operator.lshift, (v1, v2), (w1, w2))

def generic_test_numpy_rshift(v1, v2, w1, w2):
    """
    Test:
        (v1 >> v2) == (w1 >> w2)
    """
    generic_test_numpy_operator(operator.rshift, (v1, v2), (w1, w2))

def generic_test_numpy_and(v1, v2, w1, w2):
    """
    Test:
        (v1 & v2) == (w1 & w2)
    """
    generic_test_numpy_operator(operator.and, (v1, v2), (w1, w2))

def generic_test_numpy_xor(v1, v2, w1, w2):
    """
    Test:
        (v1 ^ v2) == (w1 ^ w2)
    """
    generic_test_numpy_operator(operator.xor, (v1, v2), (w1, w2))

def generic_test_numpy_or(v1, v2, w1, w2):
    """
    Test:
        (v1 | v2) == (w1 | w2)
    """
    generic_test_numpy_operator(operator.or, (v1, v2), (w1, w2))

def generic_test_numpy_(v1, v2, w1, w2):
    """
    Test:
        (v1  v2) == (w1  w2)
    """
    generic_test_numpy_operator(operator., (v1, v2), (w1, w2))

def generic_test_numpy_(v1, v2, w1, w2):
    """
    Test:
        (v1  v2) == (w1  w2)
    """
    generic_test_numpy_operator(operator., (v1, v2), (w1, w2))

def generic_test_numpy_(v1, v2, w1, w2):
    """
    Test:
        (v1  v2) == (w1  w2)
    """
    generic_test_numpy_operator(operator., (v1, v2), (w1, w2))
