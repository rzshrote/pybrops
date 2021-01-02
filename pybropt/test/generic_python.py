import inspect
import pytest

# def generic_test_function(fn, vargs, vkwargs, wargs, wkwargs):
#     pass

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
    for m in met:
        with pytest.raises(NotImplementedError):
            fn = getattr(obj, m)
            p = inspect.signature(fn).parameters
            l = list(m for m,n in p.items() if str(n)[0] != '*')
            kwargs = dict.fromkeys(l, [None] * len(l))
            fn(**kwargs)
