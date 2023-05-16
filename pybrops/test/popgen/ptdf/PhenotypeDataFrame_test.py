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

from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame
from pybrops.popgen.ptdf.PhenotypeDataFrame import is_PhenotypeDataFrame
from pybrops.popgen.ptdf.PhenotypeDataFrame import check_is_PhenotypeDataFrame


################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def df():
    yield PhenotypeDataFrame()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(PhenotypeDataFrame)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(PhenotypeDataFrame, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_df_is_abstract():
    assert_abstract_property(PhenotypeDataFrame, "col_analysis_type")

def test_ncol_is_abstract():
    assert_abstract_property(PhenotypeDataFrame, "col_analysis_effect")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_PhenotypeDataFrame_is_concrete():
    assert_concrete_function(is_PhenotypeDataFrame)

def test_check_is_PhenotypeDataFrame_is_concrete():
    assert_concrete_function(check_is_PhenotypeDataFrame)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_PhenotypeDataFrame(df):
    assert is_PhenotypeDataFrame(df)

def test_check_is_PhenotypeDataFrame(df):
    with not_raises(TypeError):
        check_is_PhenotypeDataFrame(df, "df")
    with pytest.raises(TypeError):
        check_is_PhenotypeDataFrame(None, "df")
