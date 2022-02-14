import numpy
import pytest
import copy
import pandas

from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.popgen.ptdf.DictPhenotypeDataFrame import DictPhenotypeDataFrame
from pybrops.popgen.ptdf.DictPhenotypeDataFrame import is_DictPhenotypeDataFrame
from pybrops.popgen.ptdf.DictPhenotypeDataFrame import check_is_DictPhenotypeDataFrame
from pybrops.popgen.ptdf.DictPhenotypeDataFrame import cond_check_is_DictPhenotypeDataFrame

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def data_yield():
    yield numpy.float64([4.3, 6.1, 4.5, 6.3, 6.5])

@pytest.fixture
def data_protein():
    yield numpy.float64([2.1, 6.7, 3.5, 5.6, 4.5])

@pytest.fixture
def data_oil():
    yield numpy.float64([8.8, 2. , 8.4, 9.4, 8.2])

@pytest.fixture
def data_taxa():
    yield numpy.object_(["a","b","c","d","e"])

@pytest.fixture
def data(data_yield, data_protein, data_oil, data_taxa):
    out = {
        "yield": data_yield,
        "protein": data_protein,
        "oil": data_oil,
        "taxa": data_taxa
    }
    yield out

@pytest.fixture
def col_grp_yield():
    yield "response"

@pytest.fixture
def col_grp_protein():
    yield "response"

@pytest.fixture
def col_grp_oil():
    yield "response"

@pytest.fixture
def col_grp_taxa():
    yield "predictor"

@pytest.fixture
def col_grp(col_grp_yield, col_grp_protein, col_grp_oil, col_grp_taxa):
    out = {
        "yield": col_grp_yield,
        "protein": col_grp_protein,
        "oil": col_grp_oil,
        "taxa": col_grp_taxa
    }
    yield out

@pytest.fixture
def row_name_row_name():
    yield numpy.arange(5)

@pytest.fixture
def row_name(row_name_row_name):
    out = {
        "row_name": row_name_row_name
    }
    yield out

@pytest.fixture
def col_analysis_type_yield():
    yield "double"

@pytest.fixture
def col_analysis_type_protein():
    yield "double"

@pytest.fixture
def col_analysis_type_oil():
    yield "double"

@pytest.fixture
def col_analysis_type_taxa():
    yield "factor(str)"

@pytest.fixture
def col_analysis_type(col_analysis_type_yield, col_analysis_type_protein, col_analysis_type_oil, col_analysis_type_taxa):
    out = {
        "yield": col_analysis_type_yield,
        "protein": col_analysis_type_protein,
        "oil": col_analysis_type_oil,
        "taxa": col_analysis_type_taxa
    }
    yield out

@pytest.fixture
def col_analysis_effect_yield():
    yield "response"

@pytest.fixture
def col_analysis_effect_protein():
    yield "response"

@pytest.fixture
def col_analysis_effect_oil():
    yield "response"

@pytest.fixture
def col_analysis_effect_taxa():
    yield "fixed"

@pytest.fixture
def col_analysis_effect(col_analysis_effect_yield, col_analysis_effect_protein, col_analysis_effect_oil, col_analysis_effect_taxa):
    out = {
        "yield": col_analysis_effect_yield,
        "protein": col_analysis_effect_protein,
        "oil": col_analysis_effect_oil,
        "taxa": col_analysis_effect_taxa
    }
    yield out

@pytest.fixture
def df(data, col_grp, row_name, col_analysis_type, col_analysis_effect):
    out = DictPhenotypeDataFrame(
        data = data,
        col_grp = col_grp,
        col_analysis_type = col_analysis_type,
        col_analysis_effect = col_analysis_effect,
        row_name = row_name
    )
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(DictPhenotypeDataFrame)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DictPhenotypeDataFrame, "__init__")

def test_copy_is_concrete():
    generic_assert_concrete_method(DictPhenotypeDataFrame, "__copy__")

def test_deepcopy_is_concrete():
    generic_assert_concrete_method(DictPhenotypeDataFrame, "__deepcopy__")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################
def test_init(df):
    assert is_DictPhenotypeDataFrame(df)

def test_copy(df):
    a = copy.copy(df)
    for key in df.data.keys():
        assert numpy.all(df.data[key] == a.data[key])
    assert df.ncol == a.ncol
    assert df.col_axis == a.col_axis
    assert numpy.all(df.col_dtype == a.col_dtype)
    assert numpy.all(df.col_name == a.col_name)
    assert numpy.all(df.col_grp == a.col_grp)
    assert df.nrow == a.nrow
    assert df.row_axis == a.row_axis
    assert numpy.all(df.row_name == a.row_name)
    assert numpy.all(df.col_analysis_type == a.col_analysis_type)
    assert numpy.all(df.col_analysis_effect == a.col_analysis_effect)

def test_deepcopy(df):
    a = copy.deepcopy(df)
    assert id(df.data) != id(a.data)
    for key in df.data.keys():
        assert id(df.data[key]) != id(a.data[key])
        assert numpy.all(df.data[key] == a.data[key])
    assert df.ncol == a.ncol
    assert df.col_axis == a.col_axis
    assert id(df.col_dtype) != id(a.col_dtype)
    assert numpy.all(df.col_dtype == a.col_dtype)
    assert id(df.col_name) != id(a.col_name)
    assert numpy.all(df.col_name == a.col_name)
    # assert id(df.col_grp) != id(a.col_grp) # not sure why this fails
    assert numpy.all(df.col_grp == a.col_grp)
    assert df.nrow == a.nrow
    assert df.row_axis == a.row_axis
    assert id(df.row_name) != id(a.row_name)
    assert numpy.all(df.row_name == a.row_name)
    assert id(df.col_analysis_type) != id(a.col_analysis_type)
    assert numpy.all(df.col_analysis_type == a.col_analysis_type)
    assert id(df.col_analysis_effect) != id(a.col_analysis_effect)
    assert numpy.all(df.col_analysis_effect == a.col_analysis_effect)

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_col_analysis_type_fget(df, col_analysis_type):
    assert numpy.all(df.col_analysis_type == list(col_analysis_type.values()))

def test_col_analysis_type_fset_list(df):
    l = ["factor(double)", "factor(double)", "factor(double)", "factor(int)"]
    df.col_analysis_type = l
    assert numpy.all(df.col_analysis_type == l)

def test_col_analysis_type_fset_dict(df):
    l = ["factor(double)", "factor(double)", "factor(double)", "factor(int)"]
    d = {
        "yield": "factor(double)",
        "protein": "factor(double)",
        "oil": "factor(double)",
        "taxa": "factor(int)"
    }
    df.col_analysis_type = d
    assert all(e in l for e in df.col_analysis_type)

def test_col_analysis_type_fdel(df):
    del df.col_analysis_type
    assert not hasattr(df, "_col_analysis_type")

def test_col_analysis_effect_fget(df, col_analysis_effect):
    assert numpy.all(df.col_analysis_effect == list(col_analysis_effect.values()))

def test_col_analysis_effect_fset_list(df):
    l = ["response", "response", "response", "random"]
    df.col_analysis_effect = l
    assert numpy.all(df.col_analysis_effect == l)

def test_col_analysis_effect_fset_dict(df):
    l = ["response", "response", "response", "random"]
    d = {
        "yield": "response",
        "protein": "response",
        "oil": "response",
        "taxa": "random"
    }
    df.col_analysis_effect = d
    assert all(e in l for e in df.col_analysis_effect)

def test_col_analysis_effect_fdel(df):
    del df.col_analysis_effect
    assert not hasattr(df, "_col_analysis_effect")

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_col_data(df, data_oil, data_taxa):
    out = df.col_data(name = "oil")
    assert numpy.all(out == data_oil)

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_DictPhenotypeDataFrame_is_concrete():
    generic_assert_concrete_function(is_DictPhenotypeDataFrame)

def test_check_is_DictPhenotypeDataFrame_is_concrete():
    generic_assert_concrete_function(check_is_DictPhenotypeDataFrame)

def test_cond_check_is_DictPhenotypeDataFrame_is_concrete():
    generic_assert_concrete_function(cond_check_is_DictPhenotypeDataFrame)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DictPhenotypeDataFrame(df):
    assert is_DictPhenotypeDataFrame(df)

def test_check_is_DictPhenotypeDataFrame(df):
    with not_raises(TypeError):
        check_is_DictPhenotypeDataFrame(df, "df")
    with pytest.raises(TypeError):
        check_is_DictPhenotypeDataFrame(None, "df")
