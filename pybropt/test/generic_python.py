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
    except ForbiddenException:
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

def generic_assert_docstring(obj):
    """
    Test for the presence of a docstring in an object.

    obj : Any object
        Any Python object to test for a docstring.
    """
    # test for having a docstring
    assert hasattr(obj, "__doc__")      # make sure we have the attribute
    assert isinstance(obj.__doc__, str) # make sure attribute is a string
    assert len(obj.__doc__) > 0         # make sure docstring is not empty

def generic_assert_raise_NotImplementedError(fn):
    """
    Test that a function raises NotImplementedError.

    Parameters
    ----------
    fn : function
        Function to test.
    """
    # test for raises NotImplementedError
    with pytest.raises(NotImplementedError):
        # get signature parameters
        p = inspect.signature(fn).parameters

        # get parameters as keyword arguments
        kwargs = dict((m,None) for m,n in p.items() if str(n)[0] != '*')

        # test for raises NotImplementedError
        fn(**kwargs)

def generic_assert_not_raise_NotImplementedError(fn):
    """
    Test that a function does not raise NotImplementedError.

    Parameters
    ----------
    fn : function
        Function to test.
    """
    # test for raises NotImplementedError
    with not_raises(NotImplementedError):
        # get signature parameters
        p = inspect.signature(fn).parameters

        # get parameters as keyword arguments
        kwargs = dict((m,None) for m,n in p.items() if str(n)[0] != '*')

        # test for not raise NotImplementedError
        fn(**kwargs)

def generic_assert_abstract_function(fn):
    """
    Assert an abstract function for several attributes:
        1) have a docstring
        2) raise NotImplementedError.

    Parameters
    ----------
    fn : function
        An abstract method.
    """
    generic_assert_docstring(fn)                    # assert for having a docstring
    generic_assert_raise_NotImplementedError(fn)    # assert that function is abstract

def generic_assert_concrete_function(fn):
    """
    Assert a concrete function for several attributes:
        1) have a docstring
        2) not raise NotImplementedError.

    Parameters
    ----------
    fn : function
        An abstract method.
    """
    generic_assert_docstring(fn)                        # assert for having a docstring
    generic_assert_not_raise_NotImplementedError(fn)    # assert that function is abstract

def generic_assert_hasattr(obj, a):
    """
    Assert an object has an attribute.

    Parameters
    ----------
    obj : any object
        Any Python object.
    a : str
        String of the attribute.
    """
    assert hasattr(obj, a)  # assert the object has the attribute

def generic_assert_abstract_method(obj, met):
    """
    Assert an object has an abstract method. Must have several attributes:
        1) have the method
        2) have a docstring for the method
        3) method must raise NotImplementedError

    Parameters
    ----------
    obj : any object
        Any Python object.
    met : str
        Name of the method to test
    """
    generic_assert_hasattr(obj, met)        # assert the method exists
    fn = getattr(obj, met)                  # get the method
    generic_assert_abstract_function(fn)    # assert the method is abstract

def generic_assert_concrete_method(obj, met):
    """
    Assert an object has a concrete method. Must have several attributes:
        1) have the method
        2) have a docstring for the method
        3) method must not raise NotImplementedError

    Parameters
    ----------
    obj : any object
        Any Python object.
    met : str
        Name of the method to test
    """
    generic_assert_hasattr(obj, met)        # assert the method exists
    fn = getattr(obj, met)                  # get the method
    generic_assert_concrete_function(fn)    # assert the method is abstract

def generic_assert_abstract_property(obj, prop):
    """
    Assert an object has an abstract property. Must have several attributes:
        1) have the property
        2) have a docstring for the property
        3) fget, fset, fdel methods must be abstract

    Parameters
    ----------
    obj : any object
        Any Python object.
    met : str
        Name of the method to test
    """
    generic_assert_hasattr(obj, prop)           # assert the property exists
    p = getattr(obj, prop)                      # get the property
    generic_assert_docstring(p)                 # assert the property has a docstring
    generic_assert_abstract_method(p, "fget")   # assert fget is abstract
    generic_assert_abstract_method(p, "fset")   # assert fset is abstract
    generic_assert_abstract_method(p, "fdel")   # assert fdel is abstract

def generic_test_abstract_methods(obj, mets):
    """
    Note: this is depricated.
    Test all methods for raise NotImplementedError in an object.

    obj : Any object
        Any Python object to test for method attributes.
    mets : list
        List of str for attributes to be tested.
    """
    for m in mets:
        with pytest.raises(NotImplementedError):
            # get function
            fn = getattr(obj, m)

            # get signature parameters
            p = inspect.signature(fn).parameters

            # get parameters as keyword arguments
            kwargs = dict((m,None) for m,n in p.items() if str(n)[0] != '*')

            # test for raises NotImplementedError
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
