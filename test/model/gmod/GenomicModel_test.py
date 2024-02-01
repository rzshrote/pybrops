import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.gmod.GenomicModel import check_is_GenomicModel
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmod():
    yield DummyGenomicModel()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(GenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(GenomicModel, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_model_name_is_abstract():
    assert_abstract_property(GenomicModel, "model_name")

def test_params_is_abstract():
    assert_abstract_property(GenomicModel, "hyperparams")

def test_trait_is_abstract():
    assert_abstract_property(GenomicModel, "trait")

def test_ntrait_is_abstract():
    assert_abstract_property(GenomicModel, "ntrait")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_fit_numpy_is_abstract():
    assert_abstract_method(GenomicModel, "fit_numpy")

def test_fit_is_abstract():
    assert_abstract_method(GenomicModel, "fit")

def test_predict_numpy_is_abstract():
    assert_abstract_method(GenomicModel, "predict_numpy")

def test_predict_is_abstract():
    assert_abstract_method(GenomicModel, "predict")

def test_score_numpy_is_abstract():
    assert_abstract_method(GenomicModel, "score_numpy")

def test_score_is_abstract():
    assert_abstract_method(GenomicModel, "score")

def test_gebv_numpy_is_abstract():
    assert_abstract_method(GenomicModel, "gebv_numpy")

def test_gebv_is_abstract():
    assert_abstract_method(GenomicModel, "gebv")

def test_var_G_numpy_is_abstract():
    assert_abstract_method(GenomicModel, "var_G_numpy")

def test_var_G_is_abstract():
    assert_abstract_method(GenomicModel, "var_G")

def test_var_A_numpy_is_abstract():
    assert_abstract_method(GenomicModel, "var_A_numpy")

def test_var_A_is_abstract():
    assert_abstract_method(GenomicModel, "var_A")

def test_var_a_numpy_is_abstract():
    assert_abstract_method(GenomicModel, "var_a_numpy")

def test_var_a_is_abstract():
    assert_abstract_method(GenomicModel, "var_a")

def test_bulmer_numpy_is_abstract():
    assert_abstract_method(GenomicModel, "bulmer_numpy")

def test_bulmer_is_abstract():
    assert_abstract_method(GenomicModel, "bulmer")

def test_usl_numpy_is_abstract():
    assert_abstract_method(GenomicModel, "usl_numpy")

def test_usl_is_abstract():
    assert_abstract_method(GenomicModel, "usl")

def test_lsl_numpy_is_abstract():
    assert_abstract_method(GenomicModel, "lsl_numpy")

def test_lsl_is_abstract():
    assert_abstract_method(GenomicModel, "lsl")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_GenomicModel_is_concrete():
    assert_concrete_function(check_is_GenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GenomicModel(gmod):
    with not_raises(TypeError):
        check_is_GenomicModel(gmod, "gmod")
    with pytest.raises(TypeError):
        check_is_GenomicModel(None, "gmod")
