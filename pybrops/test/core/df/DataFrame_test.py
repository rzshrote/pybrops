import inspect
import pytest

from pybrops.test.assert_python import assert_abstract_methods
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.df.DataFrame import DataFrame
from pybrops.core.df.DataFrame import check_is_DataFrame


################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def df():
    yield DataFrame()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(DataFrame)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DataFrame, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_df_is_abstract():
    assert_abstract_property(DataFrame, "data")

def test_ncol_is_abstract():
    assert_abstract_property(DataFrame, "ncol")

def test_col_axis_is_abstract():
    assert_abstract_property(DataFrame, "col_axis")

def test_col_dtype_is_abstract():
    assert_abstract_property(DataFrame, "col_dtype")

def test_col_name_is_abstract():
    assert_abstract_property(DataFrame, "col_name")

def test_col_ctype_is_abstract():
    assert_abstract_property(DataFrame, "col_grp")

def test_nrow_is_abstract():
    assert_abstract_property(DataFrame, "nrow")

def test_row_axis_is_abstract():
    assert_abstract_property(DataFrame, "row_axis")

def test_row_name_is_abstract():
    assert_abstract_property(DataFrame, "row_name")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_col_data_is_abstract():
    assert_abstract_method(DataFrame, "col_data")

def test_to_pandas_df_is_abstract():
    assert_abstract_method(DataFrame, "to_pandas_df")

def test_to_dict_is_abstract():
    assert_abstract_method(DataFrame, "to_dict")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_DataFrame_is_concrete():
    assert_concrete_function(check_is_DataFrame)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DataFrame(df):
    with not_raises(TypeError):
        check_is_DataFrame(df, "df")
    with pytest.raises(TypeError):
        check_is_DataFrame(None, "df")
