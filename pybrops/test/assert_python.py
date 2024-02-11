from abc import ABCMeta
import inspect
from typing import Callable, Generator
from contextlib import contextmanager

@contextmanager
def raises(*ExpectedExceptions: tuple[Exception]) -> Generator:
    """
    Ensure that a function raises an expected exception or one of an expected 
    set of exceptions.
    
    Parameters
    ----------
    ExpectedExceptions : tuple
        A tuple of expected exception(s).
    """
    # try to do anything in ``with`` context
    try:
        yield

    # catch any expected exceptions and ignore them as they are expected
    except ExpectedExceptions as e:
        pass

    # catch any unexpected exceptions and raise assertion error
    except Exception as e:
        raise AssertionError("unexpected exception {0} raised".format(type(e).__name__))

    # if no exceptions raised, raise an assertion error
    else:
        names = ", ".join(e.__name__ for e in ExpectedExceptions)
        raise AssertionError("did not raise any of the following exceptions: {0}".format(names))

@contextmanager
def not_raises(*ForbiddenExceptions) -> Generator:
    """
    Ensure that method does not raise an ForbiddenException.

    Parameters
    ----------
    ForbiddenException : Exception
        Forbidden Exception.
    """
    # try to do anything in ``with`` context
    try:
        yield

    # catch any forbidden exceptions and raise assertion error as they should 
    # not exist
    except ForbiddenExceptions as e:
        raise AssertionError("forbidden exception {0} raised".format(type(e).__name__))

    # catch any remaining exceptions and ignore them since we only care about 
    # whether a forbidden exception was raised
    except Exception:
        pass

def assert_docstring(obj: object) -> None:
    """
    Test for the presence of a docstring in an object.

    Parameters
    ----------
    obj : object
        Any Python object to test for a docstring.
    """
    # make sure we have the attribute
    assert hasattr(obj, "__doc__")
    
    # make sure attribute is a string
    assert isinstance(obj.__doc__, str)
    
    # make sure docstring is not empty
    assert len(obj.__doc__) > 0

def assert_function_documentation(fn: Callable) -> None:
    """
    Test a function for complete documentation.

    Parameters
    ----------
    fn : function
        A function for which to test documentation.
    """
    # assert input is a function
    assert inspect.isfunction(fn)

    # get function signature
    fn_signature = inspect.signature(fn)

    # test each parameter for type hint
    for param, hint in fn_signature.parameters.items():
        if ":" not in hint:
            raise AssertionError("parameter {0} does not have a type hint".format(param))

    # make sure we have a docstring attribute
    assert hasattr(fn, "__doc__")

    # make sure our docstring is a string
    assert isinstance(fn.__doc__, str)

    # make sure our docstring length is greater than zero
    assert len(fn.__doc__) > 0

    # make sure our parameters are in the docstring
    for param in fn_signature.parameters.keys():
        if param not in fn.__doc__:
            raise AssertionError("parameter {0} not present in docstring".format(param))

def assert_function_raises_NotImplementedError(fn: Callable) -> None:
    """
    Test that a function raises NotImplementedError.

    Parameters
    ----------
    fn : function
        Function to test.
    """
    # test for raises NotImplementedError
    with raises(NotImplementedError):
        # get signature parameters
        p = inspect.signature(fn).parameters

        # get parameters as keyword arguments
        kwargs = dict((m,None) for m,n in p.items() if str(n)[0] != '*')

        # test for raises NotImplementedError
        fn(**kwargs)

def assert_function_not_raises_NotImplementedError(fn: Callable) -> None:
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

def assert_abstract_function(fn) -> None:
    """
    Assert an abstract function for several attributes:

    1) have a docstring
    2) raise NotImplementedError.

    Parameters
    ----------
    fn : function
        An abstract method.
    """
    assert_docstring(fn)                    # assert for having a docstring
    assert_function_raises_NotImplementedError(fn)    # assert that function is abstract

def assert_concrete_function(fn) -> None:
    """
    Assert a concrete function for several attributes:

    1) have a docstring
    2) not raise NotImplementedError.

    Parameters
    ----------
    fn : function
        An abstract method.
    """
    assert_docstring(fn)                        # assert for having a docstring
    assert_function_not_raises_NotImplementedError(fn)    # assert that function is abstract

def assert_hasattr(obj: object, attr: str) -> None:
    """
    Assert an object has an attribute.

    Parameters
    ----------
    obj : object
        Any Python object.
    attr : str
        String of the attribute.
    """
    assert hasattr(obj, attr)  # assert the object has the attribute

def assert_abstract_class(obj: type) -> None:
    """
    Assert an object type is abstract. Must have several attributes:

    1) Must not have a defined ``__init__`` method within the class.
    2) Must be an ABCMeta type.
    3) Must have abstract methods.
    4) Must have a docstring for the class.

    Parameters
    ----------
    obj : type
        A Python object type.
    """
    assert '__init__' not in vars(obj)
    assert type(obj) == ABCMeta
    assert inspect.isabstract(obj)
    assert_docstring(obj)

def assert_semiabstract_class(obj: type) -> None:
    """
    Assert an object type is abstract. Must have several attributes:

    1) Must be an ABCMeta type.
    2) Must have abstract methods.
    3) Must have a docstring for the class.

    Parameters
    ----------
    obj : type
        A Python object type.
    """
    assert type(obj) == ABCMeta
    assert inspect.isabstract(obj)
    assert_docstring(obj)

def assert_mixin_class(obj: type) -> None:
    """
    Assert an object type is a mixin. Must have several attributes:

    1) Must not have a defined ``__init__`` method within the class.
    2) Must have a docstring for the class.

    Parameters
    ----------
    obj : type
        A Python object type.
    """
    assert '__init__' not in vars(obj)
    assert_docstring(obj)

def assert_concrete_class(obj: type) -> None:
    """
    Assert an object type is concrete. Must have several attributes:

    1) May inherit from ABCMeta, but cannot have abstract methods.
    2) Must have a docstring for the class.

    Parameters
    ----------
    obj : type
        A Python object type.
    """
    assert not inspect.isabstract(obj)
    assert_docstring(obj)

def assert_abstract_method(obj: object, met: str) -> None:
    """
    Assert an object has an abstract method. Must have several attributes:

    1) have the method
    2) have a docstring for the method
    3) method must raise NotImplementedError

    Parameters
    ----------
    obj : object
        Any Python object.
    met : str
        Name of the method to test
    """
    assert_hasattr(obj, met)        # assert the method exists
    fn = getattr(obj, met)                  # get the method
    assert_abstract_function(fn)    # assert the method is abstract

def assert_concrete_method(obj: object, met: str) -> None:
    """
    Assert an object has a concrete method. Must have several attributes:

    1) have the method
    2) have a docstring for the method
    3) method must not raise NotImplementedError

    Parameters
    ----------
    obj : object
        Any Python object.
    met : str
        Name of the method to test
    """
    assert_hasattr(obj, met)        # assert the method exists
    fn = getattr(obj, met)                  # get the method
    assert_concrete_function(fn)    # assert the method is abstract

def assert_abstract_property(obj: type, prop: str) -> None:
    """
    Assert an object has an abstract property. Must have several attributes:

    1) have the property
    2) have a docstring for the property
    3) fget, fset, fdel methods must be abstract if they are defined.

    Parameters
    ----------
    obj : object
        Any Python object.
    met : str
        Name of the method to test
    """
    assert_hasattr(obj, prop)               # assert the property exists
    p = getattr(obj, prop)                  # get the property
    assert_docstring(p)                     # assert the property has a docstring
    if hasattr(p, "fget") and (getattr(p, "fget") is not None):
        assert_abstract_method(p, "fget")   # assert fget is abstract
    if hasattr(p, "fset") and (getattr(p, "fset") is not None):
        assert_abstract_method(p, "fset")   # assert fset is abstract
    if hasattr(p, "fdel") and (getattr(p, "fdel") is not None):
        assert_abstract_method(p, "fdel")   # assert fdel is abstract

def assert_concrete_property_fget(obj: type, prop: str) -> None:
    """
    Assert an object has an concrete property. Must have several attributes:

    1) have the property
    2) have a docstring for the property
    3) fget, fset, fdel methods must be concrete

    Parameters
    ----------
    obj : object
        Any Python object.
    met : str
        Name of the method to test
    """
    assert_hasattr(obj, prop)           # assert the property exists
    p = getattr(obj, prop)                      # get the property
    assert_docstring(p)                 # assert the property has a docstring
    assert_concrete_method(p, "fget")   # assert fget is concrete

def assert_concrete_property(obj: type, prop: str) -> None:
    """
    Assert an object has an concrete property. Must have several attributes:

    1) have the property
    2) have a docstring for the property
    3) fget, fset, fdel methods must be concrete if they are defined.

    Parameters
    ----------
    obj : object
        Any Python object.
    met : str
        Name of the method to test
    """
    assert_hasattr(obj, prop)               # assert the property exists
    p = getattr(obj, prop)                  # get the property
    assert_docstring(p)                     # assert the property has a docstring
    if hasattr(p, "fget") and (getattr(p, "fget") is not None):
        assert_concrete_method(p, "fget")   # assert fget is concrete
    if hasattr(p, "fset") and (getattr(p, "fset") is not None):
        assert_concrete_method(p, "fset")   # assert fset is concrete
    if hasattr(p, "fdel") and (getattr(p, "fdel") is not None):
        assert_concrete_method(p, "fdel")   # assert fdel is concrete
