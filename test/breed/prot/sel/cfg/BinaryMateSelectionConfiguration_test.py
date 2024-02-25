import numpy
import pytest
from pybrops.test.assert_python import assert_method_isconcrete, assert_class_documentation, assert_class_isconcrete
from pybrops.breed.prot.sel.cfg.BinaryMateSelectionConfiguration import BinaryMateSelectionConfiguration

from .common_fixtures import *

@pytest.fixture
def selcfg(
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_pgmat,
        common_xconfig_decn_binary,
        common_xconfig_xmap,
        common_rng
    ):
    out = BinaryMateSelectionConfiguration(
        ncross = common_ncross,
        nparent = common_nparent,
        nmating = common_nmating,
        nprogeny = common_nprogeny,
        pgmat = common_pgmat,
        xconfig_decn = common_xconfig_decn_binary,
        xconfig_xmap = common_xconfig_xmap,
        rng = common_rng
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_BinaryMateSelectionConfiguration_is_concrete():
    assert_class_isconcrete(BinaryMateSelectionConfiguration)

############################## Test class docstring ############################
def test_BinaryMateSelectionConfiguration_docstring():
    assert_class_documentation(BinaryMateSelectionConfiguration)

############################# Test class properties ############################

############################## Test class methods ##############################
def test_BinaryMateSelectionConfiguration_sample_xconfig_is_concrete():
    assert_method_isconcrete(BinaryMateSelectionConfiguration, "sample_xconfig")

def test_sample_xconfig(selcfg):
    original_xconfig = selcfg.xconfig
    selcfg.sample_xconfig()
    assert numpy.any(selcfg.xconfig != original_xconfig)

def test_sample_xconfig_return_xconfig(selcfg):
    original_xconfig = selcfg.xconfig
    new_xconfig = selcfg.sample_xconfig(True)
    assert isinstance(new_xconfig, numpy.ndarray)
    assert numpy.all(selcfg.xconfig == new_xconfig)
    assert numpy.any(selcfg.xconfig != original_xconfig)
