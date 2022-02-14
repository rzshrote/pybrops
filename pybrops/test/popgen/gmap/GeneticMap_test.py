import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.popgen.gmap.GeneticMap import GeneticMap
from pybrops.popgen.gmap.GeneticMap import is_GeneticMap
from pybrops.popgen.gmap.GeneticMap import check_is_GeneticMap
from pybrops.popgen.gmap.GeneticMap import cond_check_is_GeneticMap


################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmap():
    yield GeneticMap()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(GeneticMap)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(GeneticMap, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_vrnt_chrgrp_is_abstract():
    generic_assert_abstract_property(GeneticMap, "vrnt_chrgrp")

def test_vrnt_phypos_is_abstract():
    generic_assert_abstract_property(GeneticMap, "vrnt_phypos")

def test_vrnt_genpos_is_abstract():
    generic_assert_abstract_property(GeneticMap, "vrnt_genpos")

def test_vrnt_chrgrp_name_is_abstract():
    generic_assert_abstract_property(GeneticMap, "vrnt_chrgrp_name")

def test_vrnt_chrgrp_stix_is_abstract():
    generic_assert_abstract_property(GeneticMap, "vrnt_chrgrp_stix")

def test_vrnt_chrgrp_spix_is_abstract():
    generic_assert_abstract_property(GeneticMap, "vrnt_chrgrp_spix")

def test_vrnt_chrgrp_len_is_abstract():
    generic_assert_abstract_property(GeneticMap, "vrnt_chrgrp_len")

def test_spline_is_abstract():
    generic_assert_abstract_property(GeneticMap, "spline")

def test_spline_kind_is_abstract():
    generic_assert_abstract_property(GeneticMap, "spline_kind")

def test_spline_fill_value_is_abstract():
    generic_assert_abstract_property(GeneticMap, "spline_fill_value")

################################################################################
######################## Test special abstract methods #########################
################################################################################
def test_len_is_abstract():
    generic_assert_abstract_method(GeneticMap, "__len__")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_lexsort_is_abstract():
    generic_assert_abstract_method(GeneticMap, "lexsort")

def test_reorder_is_abstract():
    generic_assert_abstract_method(GeneticMap, "reorder")

def test_sort_is_abstract():
    generic_assert_abstract_method(GeneticMap, "sort")

def test_group_is_abstract():
    generic_assert_abstract_method(GeneticMap, "group")

def test_is_grouped_is_abstract():
    generic_assert_abstract_method(GeneticMap, "is_grouped")

def test_remove_is_abstract():
    generic_assert_abstract_method(GeneticMap, "remove")

def test_select_is_abstract():
    generic_assert_abstract_method(GeneticMap, "select")

def test_prune_is_abstract():
    generic_assert_abstract_method(GeneticMap, "prune")

def test_congruence_is_abstract():
    generic_assert_abstract_method(GeneticMap, "congruence")

def test_is_congruent_is_abstract():
    generic_assert_abstract_method(GeneticMap, "is_congruent")

def test_remove_discrepancies_is_abstract():
    generic_assert_abstract_method(GeneticMap, "remove_discrepancies")

def test_build_spline_is_abstract():
    generic_assert_abstract_method(GeneticMap, "build_spline")

def test_has_spline_is_abstract():
    generic_assert_abstract_method(GeneticMap, "has_spline")

def test_interp_genpos_is_abstract():
    generic_assert_abstract_method(GeneticMap, "interp_genpos")

def test_interp_gmap_is_abstract():
    generic_assert_abstract_method(GeneticMap, "interp_gmap")

def test_gdist1g_is_abstract():
    generic_assert_abstract_method(GeneticMap, "gdist1g")

def test_gdist2g_is_abstract():
    generic_assert_abstract_method(GeneticMap, "gdist2g")

def test_gdist1p_is_abstract():
    generic_assert_abstract_method(GeneticMap, "gdist1p")

def test_gdist2p_is_abstract():
    generic_assert_abstract_method(GeneticMap, "gdist2p")

def test_to_pandas_df_is_abstract():
    generic_assert_abstract_method(GeneticMap, "to_pandas_df")

def test_to_csv_is_abstract():
    generic_assert_abstract_method(GeneticMap, "to_csv")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_GeneticMap_is_concrete():
    generic_assert_concrete_function(is_GeneticMap)

def test_check_is_GeneticMap_is_concrete():
    generic_assert_concrete_function(check_is_GeneticMap)

def test_cond_check_is_GeneticMap_is_concrete():
    generic_assert_concrete_function(cond_check_is_GeneticMap)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_GeneticMap(gmap):
    assert is_GeneticMap(gmap)

def test_check_is_GeneticMap(gmap):
    with not_raises(TypeError):
        check_is_GeneticMap(gmap, "gmap")
    with pytest.raises(TypeError):
        check_is_GeneticMap(None, "gmap")
