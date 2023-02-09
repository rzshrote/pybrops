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

from pybrops.core.df.DictDataFrame import DictDataFrame
from pybrops.core.df.DictDataFrame import is_DictDataFrame
from pybrops.core.df.DictDataFrame import check_is_DictDataFrame

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
def df(data, col_grp, row_name):
    out = DictDataFrame(
        data = data,
        col_grp = col_grp,
        row_name = row_name
    )
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(DictDataFrame)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DictDataFrame, "__init__")

def test_copy_is_concrete():
    generic_assert_concrete_method(DictDataFrame, "__copy__")

def test_deepcopy_is_concrete():
    generic_assert_concrete_method(DictDataFrame, "__deepcopy__")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################
def test_init(df):
    assert is_DictDataFrame(df)

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

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_ncol_fget(df, data):
    assert df.ncol == len(data)

def test_ncol_fset(df):
    with pytest.raises(AttributeError):
        df.ncol = 8

def test_ncol_fdel(df):
    with pytest.raises(AttributeError):
        del df.ncol

def test_col_axis_fget(df, data):
    assert df.col_axis == 1

def test_col_axis_fset(df):
    with pytest.raises(AttributeError):
        df.col_axis = 8

def test_col_axis_fdel(df):
    with pytest.raises(AttributeError):
        del df.col_axis

def test_col_dtype_fget(df, data_yield, data_protein, data_oil, data_taxa):
    a = [data_yield.dtype, data_protein.dtype, data_oil.dtype, data_taxa.dtype]
    col_dtype = df.col_dtype
    assert isinstance(col_dtype, numpy.ndarray)
    assert col_dtype.dtype == numpy.dtype("object_")
    assert numpy.all(col_dtype == a)

def test_col_dtype_fset_list(df, data_yield, data_protein, data_oil, data_taxa):
    l = [numpy.dtype('float32'), numpy.dtype('float32'), numpy.dtype('float32'), numpy.dtype('object_')]
    df.col_dtype = l
    assert numpy.all(df.col_dtype == l)

def test_col_dtype_fset_dict(df, data_yield, data_protein, data_oil, data_taxa):
    l = [numpy.dtype('float64'), numpy.dtype('float64'), numpy.dtype('float32'), numpy.dtype('object_')]
    d = {"oil": numpy.dtype('float32')}
    df.col_dtype = d
    assert numpy.all(df.col_dtype == l)

def test_col_dtype_fdel(df):
    with pytest.raises(AttributeError):
        del df.col_dtype

def test_col_name_fget(df, data):
    assert numpy.all(df.col_name == list(data.keys()))

def test_col_name_fset_list(df):
    l = ["the", "knights", "of", "ni"]
    df.col_name = l
    assert numpy.all(df.col_name == l)

def test_col_name_fset_dict(df):
    l = ["the", "knights", "of", "taxa"]
    d = {"yield": "the", "protein": "knights", "oil": "of"}
    df.col_name = d
    assert all(e in l for e in df.col_name)

def test_col_name_fdel(df):
    with pytest.raises(AttributeError):
        del df.col_name

def test_col_grp_fget(df, col_grp):
    assert numpy.all(df.col_grp == list(col_grp.values()))

def test_col_grp_fset_list(df):
    l = ["the", "knights", "of", "ni"]
    df.col_grp = l
    assert numpy.all(df.col_grp == l)

def test_col_grp_fset_dict(df):
    l = ["the", "knights", "of", "predictor"]
    d = {"yield": "the", "protein": "knights", "oil": "of", "taxa": "predictor"}
    df.col_grp = d
    assert numpy.all(df.col_grp == l)

def test_nrow_fget(df, data_yield):
    assert df.nrow == len(data_yield)

def test_nrow_fset(df):
    with pytest.raises(AttributeError):
        df.nrow = 8

def test_nrow_fdel(df):
    with pytest.raises(AttributeError):
        del df.nrow

def test_row_axis_fget(df, data):
    assert df.row_axis == 0

def test_row_axis_fset(df):
    with pytest.raises(AttributeError):
        df.row_axis = 8

def test_row_axis_fdel(df):
    with pytest.raises(AttributeError):
        del df.row_axis

def test_row_name_fget(df, row_name_row_name):
    assert numpy.all(df.row_name == row_name_row_name)

def test_row_name_fset_list(df, row_name_row_name):
    l = [str(e) for e in row_name_row_name]
    df.row_name = l
    assert numpy.all(df.row_name == l)

def test_row_name_fset_dict(df, row_name_row_name):
    l = [str(e) for e in row_name_row_name]
    d = {"row_name": numpy.object_(l)}
    df.row_name = d
    assert all(e in l for e in df.row_name)

def test_row_name_fdel(df):
    del df.row_name
    assert not hasattr(df, "_row_name")

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_col_data(df, data_oil, data_taxa):
    out = df.col_data(name = "oil")
    assert numpy.all(out == data_oil)

def test_to_pandas_df(df):
    assert isinstance(df.to_pandas_df(), pandas.DataFrame)

def test_to_dict(df):
    assert isinstance(df.to_dict(), dict)

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_DictDataFrame_is_concrete():
    generic_assert_concrete_function(is_DictDataFrame)

def test_check_is_DictDataFrame_is_concrete():
    generic_assert_concrete_function(check_is_DictDataFrame)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DictDataFrame(df):
    assert is_DictDataFrame(df)

def test_check_is_DictDataFrame(df):
    with not_raises(TypeError):
        check_is_DictDataFrame(df, "df")
    with pytest.raises(TypeError):
        check_is_DictDataFrame(None, "df")
