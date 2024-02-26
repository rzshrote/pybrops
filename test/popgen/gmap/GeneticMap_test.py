import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.popgen.gmap.GeneticMap import GeneticMap
from pybrops.popgen.gmap.GeneticMap import check_is_GeneticMap
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmap():
    yield DummyGeneticMap()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(GeneticMap)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_vrnt_chrgrp_is_abstract():
    assert_property_isabstract(GeneticMap, "vrnt_chrgrp")

def test_vrnt_phypos_is_abstract():
    assert_property_isabstract(GeneticMap, "vrnt_phypos")

def test_vrnt_genpos_is_abstract():
    assert_property_isabstract(GeneticMap, "vrnt_genpos")

def test_vrnt_chrgrp_name_is_abstract():
    assert_property_isabstract(GeneticMap, "vrnt_chrgrp_name")

def test_vrnt_chrgrp_stix_is_abstract():
    assert_property_isabstract(GeneticMap, "vrnt_chrgrp_stix")

def test_vrnt_chrgrp_spix_is_abstract():
    assert_property_isabstract(GeneticMap, "vrnt_chrgrp_spix")

def test_vrnt_chrgrp_len_is_abstract():
    assert_property_isabstract(GeneticMap, "vrnt_chrgrp_len")

def test_spline_is_abstract():
    assert_property_isabstract(GeneticMap, "spline")

def test_spline_kind_is_abstract():
    assert_property_isabstract(GeneticMap, "spline_kind")

def test_spline_fill_value_is_abstract():
    assert_property_isabstract(GeneticMap, "spline_fill_value")

################################################################################
######################## Test special abstract methods #########################
################################################################################
def test_len_is_abstract():
    assert_method_isabstract(GeneticMap, "__len__")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_lexsort_is_abstract():
    assert_method_isabstract(GeneticMap, "lexsort")

def test_reorder_is_abstract():
    assert_method_isabstract(GeneticMap, "reorder")

def test_sort_is_abstract():
    assert_method_isabstract(GeneticMap, "sort")

def test_group_is_abstract():
    assert_method_isabstract(GeneticMap, "group")

def test_is_grouped_is_abstract():
    assert_method_isabstract(GeneticMap, "is_grouped")

def test_remove_is_abstract():
    assert_method_isabstract(GeneticMap, "remove")

def test_select_is_abstract():
    assert_method_isabstract(GeneticMap, "select")

def test_congruence_is_abstract():
    assert_method_isabstract(GeneticMap, "congruence")

def test_is_congruent_is_abstract():
    assert_method_isabstract(GeneticMap, "is_congruent")

def test_remove_discrepancies_is_abstract():
    assert_method_isabstract(GeneticMap, "remove_discrepancies")

def test_build_spline_is_abstract():
    assert_method_isabstract(GeneticMap, "build_spline")

def test_has_spline_is_abstract():
    assert_method_isabstract(GeneticMap, "has_spline")

def test_interp_genpos_is_abstract():
    assert_method_isabstract(GeneticMap, "interp_genpos")

def test_interp_gmap_is_abstract():
    assert_method_isabstract(GeneticMap, "interp_gmap")

def test_gdist1g_is_abstract():
    assert_method_isabstract(GeneticMap, "gdist1g")

def test_gdist2g_is_abstract():
    assert_method_isabstract(GeneticMap, "gdist2g")

def test_gdist1p_is_abstract():
    assert_method_isabstract(GeneticMap, "gdist1p")

def test_gdist2p_is_abstract():
    assert_method_isabstract(GeneticMap, "gdist2p")

def test_to_pandas_is_abstract():
    assert_method_isabstract(GeneticMap, "to_pandas")

def test_to_csv_is_abstract():
    assert_method_isabstract(GeneticMap, "to_csv")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_GeneticMap_is_concrete():
    assert_function_isconcrete(check_is_GeneticMap)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GeneticMap(gmap):
    with not_raises(TypeError):
        check_is_GeneticMap(gmap, "gmap")
    with pytest.raises(TypeError):
        check_is_GeneticMap(None, "gmap")
