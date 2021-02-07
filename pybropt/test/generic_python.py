import inspect
import pytest
from contextlib import contextmanager

@contextmanager
def not_raises(ForbiddenException):
    """
    Ensure that method does not raise an ForbiddenException.

    ForbiddenException : Exception
        Forbidden Exception.
    """
    try:
        yield
    except ExpectedException:
        raise AssertionError("{0} raised".format(ForbiddenException.__name__))
    except Exception:
        pass

def generic_test_operator(op, v, w):
    """
    Generic test an operator.
    Tests:
        op(*v) == op(*w)

    Parameters
    ----------
    op : callable
        Operator function to call.
    """
    assert op(*v) == op(*w)

def generic_test_abstract_methods(obj, met):
    """
    Test all methods for raise NotImplementedError in an object.

    obj : Any object
        Any Python object to test for method attributes.
    met : list
        List of str for attributes to be tested.
    """
    for m in met:
        with pytest.raises(NotImplementedError):
            fn = getattr(obj, m)
            p = inspect.signature(fn).parameters
            l = list(m for m,n in p.items() if str(n)[0] != '*')
            kwargs = dict.fromkeys(l, [None] * len(l))
            fn(**kwargs)

def generic_test_concrete_methods(obj, met):
    """
    Test all methods for not raising NotImplementedError in an object.

    obj : Any object
        Any Python object to test for method attributes.
    met : list
        List of str for attributes to be tested.
    """
    for m in met:
        with not_raises(NotImplementedError):
            fn = getattr(obj, m)
            p = inspect.signature(fn).parameters
            l = list(m for m,n in p.items() if str(n)[0] != '*')
            kwargs = dict.fromkeys(l, [None] * len(l))
            fn(**kwargs)
