import inspect
import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.df import DataFrame
from pybropt.core.df import is_DataFrame
from pybropt.core.df import check_is_DataFrame
from pybropt.core.df import cond_check_is_DataFrame


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
    generic_assert_docstring(DataFrame)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DataFrame, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_df_is_abstract():
    generic_assert_abstract_property(DataFrame, "df")

def test_ncol_is_abstract():
    generic_assert_abstract_property(DataFrame, "ncol")

def test_col_axis_is_abstract():
    generic_assert_abstract_property(DataFrame, "col_axis")

def test_col_dtype_is_abstract():
    generic_assert_abstract_property(DataFrame, "col_dtype")

def test_col_name_is_abstract():
    generic_assert_abstract_property(DataFrame, "col_name")

def test_col_ctype_is_abstract():
    generic_assert_abstract_property(DataFrame, "col_ctype")

def test_nrow_is_abstract():
    generic_assert_abstract_property(DataFrame, "nrow")

def test_row_axis_is_abstract():
    generic_assert_abstract_property(DataFrame, "row_axis")

def test_row_name_is_abstract():
    generic_assert_abstract_property(DataFrame, "row_name")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_col_data_is_abstract():
    generic_assert_abstract_method(DataFrame, "col_data")

def test_to_pandas_df_is_abstract():
    generic_assert_abstract_method(DataFrame, "to_pandas_df")

def test_to_dict_is_abstract():
    generic_assert_abstract_method(DataFrame, "to_dict")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_DataFrame_is_concrete():
    generic_assert_concrete_function(is_DataFrame)

def test_check_is_DataFrame_is_concrete():
    generic_assert_concrete_function(check_is_DataFrame)

def test_cond_check_is_DataFrame_is_concrete():
    generic_assert_concrete_function(cond_check_is_DataFrame)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DataFrame(df):
    assert is_DataFrame(df)

def test_check_is_DataFrame(df):
    with not_raises(TypeError):
        check_is_DataFrame(df, "df")
    with pytest.raises(TypeError):
        check_is_DataFrame(None, "df")