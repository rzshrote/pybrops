import inspect
import pytest

def helper_test_abstract_methods(obj, met):
    for m in met:
        with pytest.raises(NotImplementedError):
            fn = getattr(obj, m)
            p = inspect.signature(fn).parameters
            kwargs = dict.fromkeys(list(p), [None] * len(p))
            fn(**kwargs)
