import operator

from . import assert_ndarray_operator

def assert_ndarray_add(v1, v2, w1, w2):
    """
    Test::

        (v1 + v2) == (w1 + w2)
    """
    assert_ndarray_operator(operator.add, (v1, v2), (w1, w2))

def assert_ndarray_sub(v1, v2, w1, w2):
    """
    Test::

        (v1 - v2) == (w1 - w2)
    """
    assert_ndarray_operator(operator.sub, (v1, v2), (w1, w2))

def assert_ndarray_mul(v1, v2, w1, w2):
    """
    Test::

        (v1 * v2) == (w1 * w2)
    """
    assert_ndarray_operator(operator.mul, (v1, v2), (w1, w2))

def assert_ndarray_matmul(v1, v2, w1, w2):
    """
    Test::

        (v1 @ v2) == (w1 @ w2)
    """
    assert_ndarray_operator(operator.matmul, (v1, v2), (w1, w2))

def assert_ndarray_truediv(v1, v2, w1, w2):
    """
    Test::

        (v1 / v2) == (w1 / w2)
    """
    assert_ndarray_operator(operator.truediv, (v1, v2), (w1, w2))

def assert_ndarray_floordiv(v1, v2, w1, w2):
    """
    Test::

        (v1 // v2) == (w1 // w2)
    """
    assert_ndarray_operator(operator.floordiv, (v1, v2), (w1, w2))

# TODO: divmod

def assert_ndarray_pow(v1, v2, w1, w2):
    """
    Test::

        (v1 ** v2) == (w1 ** w2)
    """
    assert_ndarray_operator(operator.pow, (v1, v2), (w1, w2))

def assert_ndarray_lshift(v1, v2, w1, w2):
    """
    Test::

        (v1 << v2) == (w1 << w2)
    """
    assert_ndarray_operator(operator.lshift, (v1, v2), (w1, w2))

def assert_ndarray_rshift(v1, v2, w1, w2):
    """
    Test::

        (v1 >> v2) == (w1 >> w2)
    """
    assert_ndarray_operator(operator.rshift, (v1, v2), (w1, w2))

def assert_ndarray_and(v1, v2, w1, w2):
    """
    Test::

        (v1 & v2) == (w1 & w2)
    """
    assert_ndarray_operator(operator.and_, (v1, v2), (w1, w2))

def assert_ndarray_xor(v1, v2, w1, w2):
    """
    Test::

        (v1 ^ v2) == (w1 ^ w2)
    """
    assert_ndarray_operator(operator.xor, (v1, v2), (w1, w2))

def assert_ndarray_or(v1, v2, w1, w2):
    """
    Test::

        (v1 | v2) == (w1 | w2)
    """
    assert_ndarray_operator(operator.or_, (v1, v2), (w1, w2))

def assert_ndarray_radd(v1, v2, w1, w2):
    """
    Test::

        (v2 + v1) == (w2 + w1)
    """
    assert_ndarray_operator(operator.add, (v2, v1), (w2, w1))

def assert_ndarray_rsub(v1, v2, w1, w2):
    """
    Test::

        (v2 - v1) == (w2 - w1)
    """
    assert_ndarray_operator(operator.sub, (v2, v1), (w2, w1))

def assert_ndarray_rmul(v1, v2, w1, w2):
    """
    Test::

        (v2 * v1) == (w2 * w1)
    """
    assert_ndarray_operator(operator.mul, (v2, v1), (w2, w1))

def assert_ndarray_rtruediv(v1, v2, w1, w2):
    """
    Test::

        (v2 / v1) == (w2 / w1)
    """
    assert_ndarray_operator(operator.truediv, (v2, v1), (w2, w1))

def assert_ndarray_rfloordiv(v1, v2, w1, w2):
    """
    Test::

        (v2 // v1) == (w2 // w1)
    """
    assert_ndarray_operator(operator.floordiv, (v2, v1), (w2, w1))

def assert_ndarray_rmod(v1, v2, w1, w2):
    """
    Test::

        (v2 % v1) == (w2 % w1)
    """
    assert_ndarray_operator(operator.mod, (v2, v1), (w2, w1))

# TODO: rdivmod

def assert_ndarray_rlshift(v1, v2, w1, w2):
    """
    Test::

        (v2 << v1) == (w2 << w1)
    """
    assert_ndarray_operator(operator.lshift, (v2, v1), (w2, w1))

def assert_ndarray_rrshift(v1, v2, w1, w2):
    """
    Test::

        (v2 >> v1) == (w2 >> w1)
    """
    assert_ndarray_operator(operator.rshift, (v2, v1), (w2, w1))

def assert_ndarray_rand(v1, v2, w1, w2):
    """
    Test::

        (v2 & v1) == (w2 & w1)
    """
    assert_ndarray_operator(operator.and_, (v2, v1), (w2, w1))

def assert_ndarray_rxor(v1, v2, w1, w2):
    """
    Test::

        (v2 ^ v1) == (w2 ^ w1)
    """
    assert_ndarray_operator(operator.xor, (v2, v1), (w2, w1))

def assert_ndarray_ror(v1, v2, w1, w2):
    """
    Test::

        (v2 | v1) == (w2 | w1)
    """
    assert_ndarray_operator(operator.or_, (v2, v1), (w2, w1))

def assert_ndarray_iadd(v1, v2, w1, w2):
    """
    Test::

        (v1 += v2) == (w1 += w2)
    """
    assert_ndarray_operator(operator.iadd, (v1, v2), (w1, w2))

def assert_ndarray_isub(v1, v2, w1, w2):
    """
    Test::

        (v1 -= v2) == (w1 -= w2)
    """
    assert_ndarray_operator(operator.isub, (v1, v2), (w1, w2))

def assert_ndarray_imul(v1, v2, w1, w2):
    """
    Test::

        (v1 *= v2) == (w1 *= w2)
    """
    assert_ndarray_operator(operator.imul, (v1, v2), (w1, w2))

def assert_ndarray_imatmul(v1, v2, w1, w2):
    """
    Test::

        (v1 @= v2) == (w1 @= w2)
    """
    assert_ndarray_operator(operator.imatmul, (v1, v2), (w1, w2))

def assert_ndarray_itruediv(v1, v2, w1, w2):
    """
    Test::

        (v1 /= v2) == (w1 /= w2)
    """
    assert_ndarray_operator(operator.itruediv, (v1, v2), (w1, w2))

def assert_ndarray_ifloordiv(v1, v2, w1, w2):
    """
    Test::

        (v1 //= v2) == (w1 //= w2)
    """
    assert_ndarray_operator(operator.ifloordiv, (v1, v2), (w1, w2))

def assert_ndarray_imod(v1, v2, w1, w2):
    """
    Test::

        (v1 %= v2) == (w1 %= w2)
    """
    assert_ndarray_operator(operator.imod, (v1, v2), (w1, w2))

def assert_ndarray_ipow(v1, v2, w1, w2):
    """
    Test::

        (v1 **= v2) == (w1 **= w2)
    """
    assert_ndarray_operator(operator.ipow, (v1, v2), (w1, w2))

def assert_ndarray_ilshift(v1, v2, w1, w2):
    """
    Test::

        (v1 <<= v2) == (w1 <<= w2)
    """
    assert_ndarray_operator(operator.ilshift, (v1, v2), (w1, w2))

def assert_ndarray_irshift(v1, v2, w1, w2):
    """
    Test::

        (v1 >>= v2) == (w1 >>= w2)
    """
    assert_ndarray_operator(operator.irshift, (v1, v2), (w1, w2))

def assert_ndarray_iand(v1, v2, w1, w2):
    """
    Test::

        (v1 &= v2) == (w1 &= w2)
    """
    assert_ndarray_operator(operator.iand, (v1, v2), (w1, w2))

def assert_ndarray_ixor(v1, v2, w1, w2):
    """
    Test::

        (v1 ^= v2) == (w1 ^= w2)
    """
    assert_ndarray_operator(operator.ixor, (v1, v2), (w1, w2))

def assert_ndarray_ior(v1, v2, w1, w2):
    """
    Test::

        (v1 |= v2) == (w1 |= w2)
    """
    assert_ndarray_operator(operator.ior, (v1, v2), (w1, w2))
