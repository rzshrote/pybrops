import pytest

from pybrops.test import assert_docstring
from pybrops.test import assert_concrete_method

from pybrops.opt.prob.FunctionWeight import FunctionWeight, MaximizingFunctionWeight, MinimizingFunctionWeight

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def fw_wt():
    yield 4.3

@pytest.fixture
def fw_wt_imag():
    yield complex(4.3, 5.1)

@pytest.fixture
def fw_min(fw_wt):
    yield FunctionWeight(fw_wt, "min")

@pytest.fixture
def fw_max(fw_wt):
    yield FunctionWeight(fw_wt, "max")

@pytest.fixture
def minfw_min(fw_wt):
    yield MinimizingFunctionWeight(fw_wt, "min")

@pytest.fixture
def minfw_max(fw_wt):
    yield MinimizingFunctionWeight(fw_wt, "max")

@pytest.fixture
def maxfw_min(fw_wt):
    yield MaximizingFunctionWeight(fw_wt, "min")

@pytest.fixture
def maxfw_max(fw_wt):
    yield MaximizingFunctionWeight(fw_wt, "max")

################################################################################
############################## Test class docstring ############################
################################################################################
def test_FunctionWeight_docstring():
    assert_docstring(FunctionWeight)

def test_MaximizingFunctionWeight_docstring():
    assert_docstring(MaximizingFunctionWeight)

def test_MinimizingFunctionWeight_docstring():
    assert_docstring(MinimizingFunctionWeight)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_FunctionWeight_init_is_concrete():
    assert_concrete_method(FunctionWeight, "__init__")

def test_MaximizingFunctionWeight_init_is_concrete():
    assert_concrete_method(MaximizingFunctionWeight, "__init__")

def test_MinimizingFunctionWeight_init_is_concrete():
    assert_concrete_method(MinimizingFunctionWeight, "__init__")

def test_FunctionWeight_float_is_concrete():
    assert_concrete_method(FunctionWeight, "__float__")

def test_MaximizingFunctionWeight_float_is_concrete():
    assert_concrete_method(MaximizingFunctionWeight, "__float__")

def test_MinimizingFunctionWeight_float_is_concrete():
    assert_concrete_method(MinimizingFunctionWeight, "__float__")

def test_FunctionWeight_int_is_concrete():
    assert_concrete_method(FunctionWeight, "__int__")

def test_MaximizingFunctionWeight_int_is_concrete():
    assert_concrete_method(MaximizingFunctionWeight, "__int__")

def test_MinimizingFunctionWeight_int_is_concrete():
    assert_concrete_method(MinimizingFunctionWeight, "__int__")

def test_FunctionWeight_str_is_concrete():
    assert_concrete_method(FunctionWeight, "__str__")

def test_MaximizingFunctionWeight_str_is_concrete():
    assert_concrete_method(MaximizingFunctionWeight, "__str__")

def test_MinimizingFunctionWeight_str_is_concrete():
    assert_concrete_method(MinimizingFunctionWeight, "__str__")

def test_FunctionWeight_repr_is_concrete():
    assert_concrete_method(FunctionWeight, "__repr__")

def test_MaximizingFunctionWeight_repr_is_concrete():
    assert_concrete_method(MaximizingFunctionWeight, "__repr__")

def test_MinimizingFunctionWeight_repr_is_concrete():
    assert_concrete_method(MinimizingFunctionWeight, "__repr__")

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_FunctionWeight_wt_fget(fw_wt, fw_min, fw_max):
    assert fw_wt == fw_min.wt
    assert fw_wt == fw_max.wt

def test_MinimizingFunctionWeight_wt_fget(fw_wt, minfw_min, minfw_max):
    assert fw_wt == minfw_min.wt
    assert fw_wt == minfw_max.wt

def test_MaximizingFunctionWeight_wt_fget(fw_wt, maxfw_min, maxfw_max):
    assert fw_wt == maxfw_min.wt
    assert fw_wt == maxfw_max.wt

def test_FunctionWeight_wt_fset(fw_wt, fw_wt_imag, fw_min, fw_max):
    with pytest.raises(TypeError):
        fw_min.wt = fw_wt_imag
    with pytest.raises(ValueError):
        fw_min.wt = -fw_wt
    with pytest.raises(TypeError):
        fw_max.wt = fw_wt_imag
    with pytest.raises(ValueError):
        fw_max.wt = -fw_wt

def test_MinimizingFunctionWeight_wt_fset(fw_wt, fw_wt_imag, minfw_min, minfw_max):
    with pytest.raises(TypeError):
        minfw_min.wt = fw_wt_imag
    with pytest.raises(ValueError):
        minfw_min.wt = -fw_wt
    with pytest.raises(TypeError):
        minfw_max.wt = fw_wt_imag
    with pytest.raises(ValueError):
        minfw_max.wt = -fw_wt

def test_MaximizingFunctionWeight_wt_fset(fw_wt, fw_wt_imag, maxfw_min, maxfw_max):
    with pytest.raises(TypeError):
        maxfw_min.wt = fw_wt_imag
    with pytest.raises(ValueError):
        maxfw_min.wt = -fw_wt
    with pytest.raises(TypeError):
        maxfw_max.wt = fw_wt_imag
    with pytest.raises(ValueError):
        maxfw_max.wt = -fw_wt

def test_FunctionWeight_optimization_type_fget(fw_min, fw_max):
    assert fw_min.optimization_type == "min"
    assert fw_max.optimization_type == "max"

def test_MinimizingFunctionWeight_optimization_type_fget(minfw_min, minfw_max):
    assert minfw_min.optimization_type == "min"
    assert minfw_max.optimization_type == "max"

def test_MaximizingFunctionWeight_optimization_type_fget(maxfw_min, maxfw_max):
    assert maxfw_min.optimization_type == "min"
    assert maxfw_max.optimization_type == "max"

def test_FunctionWeight_optimization_type_fset(fw_min, fw_max):
    fw_min.optimization_type = "max"
    fw_max.optimization_type = "min"
    with pytest.raises(TypeError):
        fw_min.optimization_type = None
    with pytest.raises(ValueError):
        fw_min.optimization_type = "invalid"
    with pytest.raises(TypeError):
        fw_max.optimization_type = None
    with pytest.raises(ValueError):
        fw_max.optimization_type = "invalid"

def test_MinimizingFunctionWeight_optimization_type_fset(minfw_min, minfw_max):
    minfw_min.optimization_type = "max"
    minfw_max.optimization_type = "min"
    with pytest.raises(TypeError):
        minfw_min.optimization_type = None
    with pytest.raises(ValueError):
        minfw_min.optimization_type = "invalid"
    with pytest.raises(TypeError):
        minfw_max.optimization_type = None
    with pytest.raises(ValueError):
        minfw_max.optimization_type = "invalid"

def test_MaximizingFunctionWeight_optimization_type_fset(maxfw_min, maxfw_max):
    maxfw_min.optimization_type = "max"
    maxfw_max.optimization_type = "min"
    with pytest.raises(TypeError):
        maxfw_min.optimization_type = None
    with pytest.raises(ValueError):
        maxfw_min.optimization_type = "invalid"
    with pytest.raises(TypeError):
        maxfw_max.optimization_type = None
    with pytest.raises(ValueError):
        maxfw_max.optimization_type = "invalid"

################################################################################
######################### Test special class functions #########################
################################################################################
def test_FunctionWeight_float(fw_wt, fw_min, fw_max):
    assert float(fw_min) == float(fw_wt)
    assert float(fw_max) == float(fw_wt)

def test_MinimizingFunctionWeight_float(fw_wt, minfw_min, minfw_max):
    assert float(minfw_min) == float(fw_wt)
    assert float(minfw_max) == -float(fw_wt)

def test_MaximizingFunctionWeight_float(fw_wt, maxfw_min, maxfw_max):
    assert float(maxfw_min) == -float(fw_wt)
    assert float(maxfw_max) == float(fw_wt)

def test_FunctionWeight_int(fw_wt, fw_min, fw_max):
    assert int(fw_min) == int(fw_wt)
    assert int(fw_max) == int(fw_wt)

def test_MinimizingFunctionWeight_int(fw_wt, minfw_min, minfw_max):
    assert int(minfw_min) == int(fw_wt)
    assert int(minfw_max) == -int(fw_wt)

def test_MaximizingFunctionWeight_int(fw_wt, maxfw_min, maxfw_max):
    assert int(maxfw_min) == -int(fw_wt)
    assert int(maxfw_max) == int(fw_wt)

def test_FunctionWeight_str(fw_wt, fw_min, fw_max):
    assert str(fw_min) == str(fw_wt)
    assert str(fw_max) == str(fw_wt)

def test_MinimizingFunctionWeight_str(fw_wt, minfw_min, minfw_max):
    assert str(minfw_min) == str(fw_wt)
    assert str(minfw_max) == str(-fw_wt)

def test_MaximizingFunctionWeight_str(fw_wt, maxfw_min, maxfw_max):
    assert str(maxfw_min) == str(-fw_wt)
    assert str(maxfw_max) == str(fw_wt)

def test_FunctionWeight_repr(fw_wt, fw_min, fw_max):
    assert repr(fw_min) == repr(fw_wt)
    assert repr(fw_max) == repr(fw_wt)

def test_MinimizingFunctionWeight_repr(fw_wt, minfw_min, minfw_max):
    assert repr(minfw_min) == repr(fw_wt)
    assert repr(minfw_max) == repr(-fw_wt)

def test_MaximizingFunctionWeight_repr(fw_wt, maxfw_min, maxfw_max):
    assert repr(maxfw_min) == repr(-fw_wt)
    assert repr(maxfw_max) == repr(fw_wt)
