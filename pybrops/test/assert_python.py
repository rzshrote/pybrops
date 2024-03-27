from abc import ABCMeta
import inspect
import re
from typing import Callable, Generator, Tuple
from contextlib import contextmanager

@contextmanager
def raises(*ExpectedExceptions: Tuple[Exception,...]) -> Generator:
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
        raise AssertionError("unexpected {0} exception raised".format(type(e).__name__))

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
        raise AssertionError("forbidden {0} exception raised".format(type(e).__name__))

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

def assert_hasattr(obj: object, name: str) -> None:
    """
    Assert an object has an attribute.

    Parameters
    ----------
    obj : object
        Any Python object.
    attr : str
        String of the attribute.
    """
    assert hasattr(obj, name)  # assert the object has the attribute

######################### Function assertion functions #########################

def assert_function_documentation(fn: Callable) -> None:
    """
    Test a function for complete documentation.

    Parameters
    ----------
    fn : function
        A function for which to test documentation.
    """
    # assert input is callable or a function
    # this will fail if input function is a builtin function
    if not (callable(fn) and inspect.isfunction(fn)):
        raise AssertionError("object ``{0}`` is not callable and a function".format(fn.__name__))

    # make sure we have a docstring attribute
    if not hasattr(fn, "__doc__"):
        raise AssertionError("in function ``{0}``: docstring is not present".format(fn.__name__))

    # make sure our docstring is a string
    if not isinstance(fn.__doc__, str):
        raise AssertionError("in function ``{0}``: docstring is not a string".format(fn.__name__))

    # make sure our docstring length is greater than zero
    if len(fn.__doc__) <= 0:
        raise AssertionError("in function ``{0}``: docstring is empty".format(fn.__name__))

    # get function signature
    fn_signature = inspect.signature(fn)

    # get method parameters
    fn_parameters = {param:hint for param,hint in fn_signature.parameters.items() if param not in ("self","cls")}

    # get method return
    fn_return = fn_signature.return_annotation

    # test each parameter for type hint
    for param, hint in fn_parameters.items():
        if hint.annotation is hint.empty:
            raise AssertionError("in function ``{0}``: parameter ``{1}`` does not have a type hint".format(fn.__name__,param))

    # test the return type hint
    if fn_return is fn_signature.empty:
        raise AssertionError("in function ``{0}``: return type hint not present".format(fn.__name__))

    # make sure our docstring has a Parameters section
    if len(fn_parameters) > 0:
        if not re.search(r'Parameters\n\s*----------', fn.__doc__):
            raise AssertionError("in function ``{0}``: no ``Parameters`` section in docstring".format(fn.__name__))
    
    # make sure our docstring has a Returns section
    if fn_return is not None:
        if not re.search(r'Returns\n\s*-------', fn.__doc__):
            raise AssertionError("in function ``{0}``: no ``Returns`` section in docstring".format(fn.__name__))

    # make sure our parameters are in the docstring
    for param in fn_parameters.keys():
        if not re.search(param + r'\s*:\s*\w+[^\n]*\n', fn.__doc__):
            raise AssertionError("in function ``{0}``: parameter ``{1}`` not present in docstring".format(fn.__name__,param))

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
        params = inspect.signature(fn).parameters

        # get parameters as keyword arguments
        kwargs = dict((m,None) for m,n in params.items() if str(n)[0] != '*')

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
        params = inspect.signature(fn).parameters

        # get parameters as keyword arguments
        kwargs = dict((m,None) for m,n in params.items() if str(n)[0] != '*')

        # test for not raise NotImplementedError
        fn(**kwargs)

def assert_function_isabstract(fn: Callable) -> None:
    """
    Assert an abstract function for several attributes:

    1) have documentation
    2) raise NotImplementedError.

    Parameters
    ----------
    fn : function
        An abstract method.
    """
    # assert for having a docstring
    assert_function_documentation(fn)

    # assert that function is abstract
    assert_function_raises_NotImplementedError(fn)

def assert_function_isconcrete(fn: Callable) -> None:
    """
    Assert a concrete function for several attributes:

    1) have documentation
    2) not raise NotImplementedError.

    Parameters
    ----------
    fn : function
        An abstract method.
    """
    # assert for having a docstring
    assert_function_documentation(fn)

    # assert that function is abstract
    assert_function_not_raises_NotImplementedError(fn)

########################## Method assertion functions ##########################

def assert_method_documentation(obj: type, name: str) -> None:
    """
    Test a method for complete documentation.

    Parameters
    ----------
    obj : type
        Object type for which the method is attached.
    name : str
        Name of the method to test.
    """
    # assert inputs are correct
    if not isinstance(obj, type):
        raise AssertionError("input ``obj`` is not a type")
    
    if not isinstance(name, str):
        raise AssertionError("input ``name`` is not a string")

    # assert that the type has the method
    if not hasattr(obj, name):
        raise AssertionError("type ``{0}`` has no attribute ``{1}``".format(obj.__name__,name))

    # get the method
    met = getattr(obj, name)

    # assert that the object is a method
    # this will fail if input function is a builtin method
    if not (callable(met) and inspect.isfunction(met)):
        raise AssertionError("object ``{0}.{1}`` is not callable and a method".format(obj.__name__,name))

    # make sure we have a docstring attribute
    if not hasattr(met, "__doc__"):
        raise AssertionError("in method ``{0}.{1}``: docstring is not present".format(obj.__name__,met.__name__))

    # make sure our docstring is a string
    if not isinstance(met.__doc__, str):
        raise AssertionError("in method ``{0}.{1}``: docstring is not a string".format(obj.__name__,met.__name__))

    # make sure our docstring length is greater than zero
    if len(met.__doc__) <= 0:
        raise AssertionError("in method ``{0}.{1}``: docstring is empty".format(obj.__name__,met.__name__))

    # get method signature
    met_signature = inspect.signature(met)

    # get method parameters
    met_parameters = {param:hint for param,hint in met_signature.parameters.items() if param not in ("self","cls")}

    # get method return
    met_return = met_signature.return_annotation

    # test each parameter for type hint
    for param, hint in met_parameters.items():
        if hint.annotation is hint.empty:
            raise AssertionError("in method ``{0}.{1}``: parameter ``{2}`` does not have a type hint".format(obj.__name__,met.__name__,param))

    # test the return type hint
    if met_return is met_signature.empty:
        raise AssertionError("in method ``{0}.{1}``: return type hint not present".format(obj.__name__,met.__name__))

    # make sure our docstring has a Parameters section
    if len(met_parameters) > 0:
        if not re.search(r'Parameters\n\s*----------', met.__doc__):
            raise AssertionError("in method ``{0}.{1}``: no ``Parameters`` section in docstring".format(obj.__name__,met.__name__))
    
    # make sure our docstring has a Returns section
    if met_return is not None:
        if not re.search(r'Returns\n\s*-------', met.__doc__):
            raise AssertionError("in method ``{0}.{1}``: no ``Returns`` section in docstring".format(obj.__name__,met.__name__))

    # make sure our parameters are in the docstring
    for param in met_parameters.keys():
        if not re.search(param + r'\s*:\s*\w+[^\n]*\n', met.__doc__):
            raise AssertionError("in method ``{0}.{1}``: parameter ``{2}`` not present in docstring".format(obj.__name__,met.__name__,param))

def assert_method_raises_NotImplementedError(obj: type, name: str) -> None:
    """
    Test that a method raises NotImplementedError.

    Parameters
    ----------
    obj : type
        Object type for which the method is attached.
    name : str
        Name of the method to test.
    """
    # get the method
    met = getattr(obj, name)

    # test for raises NotImplementedError
    with raises(NotImplementedError):
        # get signature parameters
        params = inspect.signature(met).parameters

        # get parameters as keyword arguments
        kwargs = dict((m,None) for m,n in params.items() if str(n)[0] != '*')

        # test for raises NotImplementedError
        met(**kwargs)

def assert_method_not_raises_NotImplementedError(obj: type, name: str) -> None:
    """
    Test that a method does not raise NotImplementedError.

    Parameters
    ----------
    obj : type
        Object type for which the method is attached.
    name : str
        Name of the method to test.
    """
    # get the method
    met = getattr(obj, name)

    # test for raises NotImplementedError
    with not_raises(NotImplementedError):
        # get signature parameters
        params = inspect.signature(met).parameters

        # get parameters as keyword arguments
        kwargs = dict((m,None) for m,n in params.items() if str(n)[0] != '*')

        # test for not raise NotImplementedError
        met(**kwargs)

def assert_method_isabstract(obj: object, name: str) -> None:
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
    # assert the method has documentation
    assert_method_documentation(obj, name)

    # assert the method is abstract
    assert_method_raises_NotImplementedError(obj, name)

def assert_method_isconcrete(obj: object, name: str) -> None:
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
    # assert the method exists
    assert hasattr(obj, name)

    # assert the method has documentation
    assert_method_documentation(obj, name)

    # assert the method is concrete
    assert_method_not_raises_NotImplementedError(obj, name)

####################### Classmethod assertion functions ########################

def assert_classmethod_documentation(obj: type, name: str) -> None:
    """
    Test a classmethod for complete documentation.

    Parameters
    ----------
    obj : type
        Object type for which the classmethod is attached.
    name : str
        Name of the classmethod to test.
    """
    # assert inputs are correct
    if not isinstance(obj, type):
        raise AssertionError("input ``obj`` is not a type")
    
    if not isinstance(name, str):
        raise AssertionError("input ``name`` is not a string")

    # assert that the type has the classmethod
    if not hasattr(obj, name):
        raise AssertionError("object ``obj`` has no attribute ``{1}``".format(obj.__name__,name))

    # get the classmethod
    met = getattr(obj, name)

    # assert that the object is a method
    # this will fail if input function is a builtin method
    if not (callable(met) and inspect.ismethod(met)):
        raise AssertionError("object ``{0}.{1}`` is not callable and a classmethod".format(obj.__name__,met.__name__))

    # make sure we have a docstring attribute
    if not hasattr(met, "__doc__"):
        raise AssertionError("in classmethod ``{0}.{1}``: docstring is not present".format(obj.__name__,met.__name__))

    # make sure our docstring is a string
    if not isinstance(met.__doc__, str):
        raise AssertionError("in classmethod ``{0}.{1}``: docstring is not a string".format(obj.__name__,met.__name__))

    # make sure our docstring length is greater than zero
    if len(met.__doc__) <= 0:
        raise AssertionError("in classmethod ``{0}.{1}``: docstring is empty".format(obj.__name__,met.__name__))

    # get classmethod signature
    met_signature = inspect.signature(met)

    # get classmethod parameters
    met_parameters = {param:hint for param,hint in met_signature.parameters.items() if param not in ("self","cls")}

    # get classmethod return
    met_return = met_signature.return_annotation

    # test each parameter for type hint
    for param, hint in met_parameters.items():
        if hint.annotation is hint.empty:
            raise AssertionError("in classmethod ``{0}.{1}``: parameter ``{2}`` does not have a type hint".format(obj.__name__,met.__name__,param))

    # test the return type hint
    if met_return is met_signature.empty:
        raise AssertionError("in classmethod ``{0}.{1}``: return type hint not present".format(obj.__name__,met.__name__))

    # make sure our docstring has a Parameters section
    if len(met_parameters) > 0:
        if not re.search(r'Parameters\n\s*----------', met.__doc__):
            raise AssertionError("in classmethod ``{0}.{1}``: no ``Parameters`` section in docstring".format(obj.__name__,met.__name__))
    
    # make sure our docstring has a Returns section
    if met_return is not None:
        if not re.search(r'Returns\n\s*-------', met.__doc__):
            raise AssertionError("in classmethod ``{0}.{1}``: no ``Returns`` section in docstring".format(obj.__name__,met.__name__))

    # make sure our parameters are in the docstring
    for param in met_parameters.keys():
        if not re.search(param + r'\s*:\s*\w+[^\n]*\n', met.__doc__):
            raise AssertionError("in classmethod ``{0}.{1}``: parameter ``{2}`` not present in docstring".format(obj.__name__,met.__name__,param))

def assert_classmethod_raises_NotImplementedError(obj: type, name: str) -> None:
    """
    Test that a classmethod raises NotImplementedError.

    Parameters
    ----------
    obj : type
        Object type for which the classmethod is attached.
    name : str
        Name of the classmethod to test.
    """
    # get the classmethod
    met = getattr(obj, name)

    # test for raises NotImplementedError
    with raises(NotImplementedError):
        # get signature parameters
        params = inspect.signature(met).parameters

        # get parameters as keyword arguments
        kwargs = dict((m,None) for m,n in params.items() if str(n)[0] != '*')

        # test for raises NotImplementedError
        met(**kwargs)

def assert_classmethod_not_raises_NotImplementedError(obj: type, name: str) -> None:
    """
    Test that a classmethod does not raise NotImplementedError.

    Parameters
    ----------
    obj : type
        Object type for which the classmethod is attached.
    name : str
        Name of the classmethod to test.
    """
    # get the classmethod
    met = getattr(obj, name)

    # test for raises NotImplementedError
    with not_raises(NotImplementedError):
        # get signature parameters
        params = inspect.signature(met).parameters

        # get parameters as keyword arguments
        kwargs = dict((m,None) for m,n in params.items() if str(n)[0] != '*')

        # test for not raise NotImplementedError
        met(**kwargs)

def assert_classmethod_isabstract(obj: object, name: str) -> None:
    """
    Assert an object has an abstract classmethod. Must have several attributes:

    1) have the classmethod
    2) have a docstring for the classmethod
    3) classmethod must raise NotImplementedError

    Parameters
    ----------
    obj : object
        Any Python object.
    met : str
        Name of the classmethod to test
    """
    # assert the classmethod has documentation
    assert_classmethod_documentation(obj, name)

    # assert the classmethod is abstract
    assert_classmethod_raises_NotImplementedError(obj, name)

def assert_classmethod_isconcrete(obj: object, name: str) -> None:
    """
    Assert an object has a concrete classmethod. Must have several attributes:

    1) have the classmethod
    2) have a docstring for the classmethod
    3) classmethod must not raise NotImplementedError

    Parameters
    ----------
    obj : object
        Any Python object.
    met : str
        Name of the classmethod to test
    """
    # assert the classmethod has documentation
    assert_classmethod_documentation(obj, name)

    # assert the classmethod is concrete
    assert_classmethod_not_raises_NotImplementedError(obj, name)

######################### Property assertion functions #########################

def assert_property_documentation(obj: type, name: str) -> None:
    """
    Test a property for complete documentation.

    Parameters
    ----------
    obj : type
        Object type for which the method is attached.
    name : str
        Name of the method to test.
    """
    # assert inputs are correct
    if not isinstance(obj, type):
        raise AssertionError("input ``obj`` is not a type")
    
    if not isinstance(name, str):
        raise AssertionError("input ``name`` is not a string")

    # assert that the type has the method
    if not hasattr(obj, name):
        raise AssertionError("object ``obj`` has no attribute ``{1}``".format(obj.__name__,name))

    # get the property
    prop = getattr(obj, name)

    # assert that the object is a property
    if not isinstance(prop, property):
        raise AttributeError("object ``{0}.{1}`` is not a property".format(obj.__name__,name))

    # make sure we have a docstring attribute
    if not hasattr(prop, "__doc__"):
        raise AssertionError("in property ``{0}.{1}``: docstring is not present".format(obj.__name__,met.__name__))

    # make sure our docstring is a string
    if not isinstance(prop.__doc__, str):
        raise AssertionError("in property ``{0}.{1}``: docstring is not a string".format(obj.__name__,met.__name__))

    # make sure our docstring length is greater than zero
    if len(prop.__doc__) <= 0:
        raise AssertionError("in property ``{0}.{1}``: docstring is empty".format(obj.__name__,met.__name__))

    # test documentation for fget, fset, fdel
    for propmet in ("fget","fset","fdel"):
        # test the method of the property
        if hasattr(prop, propmet) and (getattr(prop, propmet) is not None):
            # get the method
            met = getattr(prop, propmet)

            # make sure we have a docstring attribute
            if not hasattr(met, "__doc__"):
                raise AssertionError("in property method ``{0}.{1}.{2}``: docstring is not present".format(obj.__name__,name,propmet))

            # make sure our docstring is a string
            if not isinstance(met.__doc__, str):
                raise AssertionError("in property method ``{0}.{1}.{2}``: docstring is not a string".format(obj.__name__,name,propmet))

            # make sure our docstring length is greater than zero
            if len(met.__doc__) <= 0:
                raise AssertionError("in property method ``{0}.{1}.{2}``: docstring is empty".format(obj.__name__,name,propmet))

            # assert that the object is a method
            # this will fail if input function is a builtin method
            if not (callable(met) and inspect.isfunction(met)):
                raise AssertionError("property method ``{0}.{1}.{2}`` is not callable and a function".format(obj.__name__,name,propmet))

            # get method signature
            met_signature = inspect.signature(met)

            # get method parameters
            met_parameters = {param:hint for param,hint in met_signature.parameters.items() if param not in ("self","cls")}

            # get method return
            met_return = met_signature.return_annotation

            # don't test for parameter names in docstring since these docs are 
            # typically smaller
            # 
            # make sure our docstring has a Parameters section
            # if len(met_parameters) > 0:
            #     if not re.search(r'Parameters\n\s*----------', met.__doc__):
            #         raise AssertionError("in method ``{0}.{1}``: no ``Parameters`` section in docstring".format(obj.__name__,met.__name__))
            # 
            # make sure our docstring has a Returns section
            # if met_return is not None:
            #     if not re.search(r'Returns\n\s*-------', met.__doc__):
            #         raise AssertionError("in method ``{0}.{1}``: no ``Returns`` section in docstring".format(obj.__name__,met.__name__))

            # test each parameter for type hint
            for param, hint in met_parameters.items():
                if hint.annotation is hint.empty:
                    raise AssertionError("in property ``{0}.{1}.{2}``: parameter ``{3}`` does not have a type hint".format(obj.__name__,name,propmet,param))

            # test the return type hint
            if met_return is met_signature.empty:
                raise AssertionError("in property ``{0}.{1}.{2}``: return type hint not present".format(obj.__name__,name,propmet))


def assert_property_raises_NotImplementedError(obj: type, name: str) -> None:
    """
    Test that a property raises NotImplementedError.

    Parameters
    ----------
    obj : type
        Object type for which the method is attached.
    name : str
        Name of the method to test.
    """
    # get the property
    prop = getattr(obj, name)

    # assert that the object is a property
    if not isinstance(prop, property):
        raise AttributeError("object ``{0}.{1}`` is not a property".format(obj.__name__,name))

    # test methods fget, fset, fdel
    for propmet in ("fget","fset","fdel"):
        # test the method of the property
        if hasattr(prop, propmet) and (getattr(prop, propmet) is not None):

            # get the method
            met = getattr(prop, propmet)

            # assert that the object is a method
            # this will fail if input function is a builtin method
            if not (callable(met) and inspect.isfunction(met)):
                raise AssertionError("property method ``{0}.{1}.{2}`` is not callable and a function".format(obj.__name__,name,propmet))

            # test for raises NotImplementedError
            with raises(NotImplementedError):
                # get signature parameters
                params = inspect.signature(met).parameters

                # get parameters as keyword arguments
                kwargs = dict((m,None) for m,n in params.items() if str(n)[0] != '*')

                # test for raises NotImplementedError
                met(**kwargs)

def assert_property_not_raises_NotImplementedError(obj: type, name: str) -> None:
    """
    Test that a method does not raise NotImplementedError.

    Parameters
    ----------
    obj : type
        Object type for which the method is attached.
    name : str
        Name of the method to test.
    """
    # get the property
    prop = getattr(obj, name)

    # assert that the object is a property
    if not isinstance(prop, property):
        raise AttributeError("object ``{0}.{1}`` is not a property".format(obj.__name__,name))

    # test methods fget, fset, fdel
    for propmet in ("fget","fset","fdel"):
        # test the method of the property
        if hasattr(prop, propmet) and (getattr(prop, propmet) is not None):

            # get the method
            met = getattr(prop, "fget")

            # assert that the object is a method
            # this will fail if input function is a builtin method
            if not (callable(met) and inspect.isfunction(met)):
                raise AssertionError("property method ``{0}.{1}.{2}`` is not callable and a function".format(obj.__name__,name,propmet))

            # test for raises NotImplementedError
            with not_raises(NotImplementedError):
                # get signature parameters
                params = inspect.signature(met).parameters

                # get parameters as keyword arguments
                kwargs = dict((m,None) for m,n in params.items() if str(n)[0] != '*')

                # test for raises NotImplementedError
                met(**kwargs)

def assert_property_isabstract(obj: type, name: str) -> None:
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
    # assert the property has documentation
    assert_property_documentation(obj, name)
    
    # assert the property is abstract
    assert_property_raises_NotImplementedError(obj, name)

def assert_property_isconcrete(obj: type, name: str) -> None:
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
    # assert the property has documentation
    assert_property_documentation(obj, name)
    
    # assert the property is abstract
    assert_property_not_raises_NotImplementedError(obj, name)

########################### Class assertion functions ##########################

def assert_class_documentation(obj: type) -> None:
    """
    Test a class for complete documentation.

    Parameters
    ----------
    obj : type
        Object type to test.
    """
    # assert inputs are correct
    if not isinstance(obj, type):
        raise AssertionError("input ``obj`` is not a type")

    # assert that the object is a class
    if not inspect.isclass(obj):
        raise AssertionError("type ``{0}`` is not a class".format(obj.__name__))

    # make sure we have a docstring attribute
    if not hasattr(obj, "__doc__"):
        raise AssertionError("in class ``{0}``: docstring is not present".format(obj.__name__))

    # make sure our docstring is a string
    if not isinstance(obj.__doc__, str):
        raise AssertionError("in class ``{0}``: docstring is not a string".format(obj.__name__))

    # make sure our docstring length is greater than zero
    if len(obj.__doc__) <= 0:
        raise AssertionError("in class ``{0}``: docstring is empty".format(obj.__name__))

def assert_class_isabstract(obj: type) -> None:
    """
    Assert an object type is abstract. Must have several attributes:

    1) Must not have a defined ``__init__`` method within the class.
    2) Must be an ABCMeta type.
    3) Must have abstract methods.
    4) Must have documentation for the class.

    Parameters
    ----------
    obj : type
        A Python object type.
    """
    # assert class has documentation
    assert_class_documentation(obj)
    
    # assert the abstract class is abstract
    if not inspect.isabstract(obj):
        raise AssertionError("class ``{0}`` is not abstract".format(obj.__name__))

    # assert the class does not have an __init__ function defined
    if '__init__' in vars(obj):
        raise AssertionError("in abstract class ``{0}``: the method ``__init__`` is defined".format(obj.__name__))
    
    # assert the abstract class type is ABCMeta
    if type(obj) != ABCMeta:
        raise AssertionError("in abstract class ``{0}``: type is not ``{1}``".format(obj.__name__,ABCMeta.__name__))

def assert_class_issemiabstract(obj: type) -> None:
    """
    Assert an object type is abstract. Must have several attributes:

    1) Must be an ABCMeta type.
    2) Must have abstract methods.
    3) Must have documentation for the class.

    Parameters
    ----------
    obj : type
        A Python object type.
    """
    # assert class has documentation
    assert_class_documentation(obj)
    
    # assert the abstract class is abstract
    if not inspect.isabstract(obj):
        raise AssertionError("class ``{0}`` is not semi-abstract".format(obj.__name__))

    # assert the abstract class type is ABCMeta
    if type(obj) != ABCMeta:
        raise AssertionError("in semi-abstract class ``{0}``: type is not ``{1}``".format(obj.__name__,ABCMeta.__name__))

def assert_class_ismixin(obj: type) -> None:
    """
    Assert an object type is a mixin. Must have several attributes:

    1) Must not have a defined ``__init__`` method within the class.
    2) Must have documentation for the class.

    Parameters
    ----------
    obj : type
        A Python object type.
    """
    # assert class has documentation
    assert_class_documentation(obj)
    
    # assert the class does not have an __init__ function defined
    if '__init__' in vars(obj):
        raise AssertionError("in mixin class ``{0}``: the method ``__init__`` is defined".format(obj.__name__))

def assert_class_isconcrete(obj: type) -> None:
    """
    Assert an object type is concrete. Must have several attributes:

    1) May inherit from ABCMeta, but cannot have abstract methods.
    2) Must have documentation for the class.

    Parameters
    ----------
    obj : type
        A Python object type.
    """
    # assert class has documentation
    assert_class_documentation(obj)
    
    # assert the abstract class is abstract
    if inspect.isabstract(obj):
        raise AssertionError("class ``{0}`` is abstract".format(obj.__name__))

########################## Module assertion functions ##########################

def assert_module_documentation(obj: object) -> None:
    """
    Test a module for complete documentation.

    Parameters
    ----------
    obj : object
    """
    # assert that the object is a class
    if not inspect.ismodule(obj):
        raise AssertionError("type ``{0}`` is not a module".format(obj.__name__))

    # make sure we have a docstring attribute
    if not hasattr(obj, "__doc__") or obj.__doc__ is None:
        raise AssertionError("in module ``{0}``: docstring is not present".format(obj.__name__))

    # make sure our docstring is a string
    if not isinstance(obj.__doc__, str):
        raise AssertionError("in module ``{0}``: docstring is not a string".format(obj.__name__))

    # make sure our docstring length is greater than zero
    if len(obj.__doc__) <= 0:
        raise AssertionError("in module ``{0}``: docstring is empty".format(obj.__name__))
