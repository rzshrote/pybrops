import numpy
import pytest
from pybrops.test.assert_python import assert_method_isconcrete, assert_class_documentation, assert_class_isconcrete, assert_module_documentation
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

############################ Test module attributes ############################
def test_BinaryMateSelectionConfiguration_module_is_concrete():
    import pybrops.breed.prot.sel.cfg.BinaryMateSelectionConfiguration
    assert_module_documentation(pybrops.breed.prot.sel.cfg.BinaryMateSelectionConfiguration)

############################ Test class attributes #############################
def test_BinaryMateSelectionConfiguration_is_concrete():
    assert_class_isconcrete(BinaryMateSelectionConfiguration)

##################### Test class special concrete methods ######################
def test_BinaryMateSelectionConfiguration___init___is_concrete():
    assert_method_isconcrete(BinaryMateSelectionConfiguration, "__init__")

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
