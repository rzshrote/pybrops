import numpy
import pytest
from pybrops.test.assert_python import assert_method_isconcrete, assert_module_documentation, assert_property_isconcrete, assert_class_documentation, assert_class_isconcrete, not_raises
from pybrops.breed.prot.sel.cfg.BinarySelectionConfiguration import BinarySelectionConfiguration

from .common_fixtures import *

@pytest.fixture
def selcfg(
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_pgmat,
        common_xconfig_decn_binary,
        common_rng
    ):
    out = BinarySelectionConfiguration(
        ncross = common_ncross,
        nparent = common_nparent,
        nmating = common_nmating,
        nprogeny = common_nprogeny,
        pgmat = common_pgmat,
        xconfig_decn = common_xconfig_decn_binary,
        rng = common_rng
    )
    yield out

############################ Test module attributes ############################
def test_BinarySelectionConfiguration_module_is_concrete():
    import pybrops.breed.prot.sel.cfg.BinarySelectionConfiguration
    assert_module_documentation(pybrops.breed.prot.sel.cfg.BinarySelectionConfiguration)

############################ Test class attributes #############################
def test_BinarySelectionConfiguration_is_concrete():
    assert_class_isconcrete(BinarySelectionConfiguration)

##################### Test class special concrete methods ######################
def test_BinarySelectionConfiguration___init___is_concrete():
    assert_method_isconcrete(BinarySelectionConfiguration, "__init__")

############################# Test class properties ############################

### xconfig_decn ###
def test_BinarySelectionConfiguration_xconfig_decn_is_concrete():
    assert_property_isconcrete(BinarySelectionConfiguration, "xconfig_decn")

def test_xconfig_decn_fget(selcfg, common_xconfig_decn_binary):
    assert numpy.all(selcfg.xconfig_decn == common_xconfig_decn_binary)

def test_xconfig_decn_fset(selcfg, common_xconfig_decn_bool, common_xconfig_decn_binary, common_nconfig):
    with not_raises(Exception):
        selcfg.xconfig_decn = common_xconfig_decn_bool
    with not_raises(Exception):
        selcfg.xconfig_decn = common_xconfig_decn_binary
    with not_raises(Exception):
        selcfg.xconfig_decn = numpy.random.randint(0, 2, common_nconfig)

def test_xconfig_decn_fset_TypeError(selcfg, common_xconfig_decn_real, common_nconfig):
    with pytest.raises(TypeError):
        selcfg.xconfig_decn = object()
    with pytest.raises(TypeError):
        selcfg.xconfig_decn = None
    with pytest.raises(TypeError):
        selcfg.xconfig_decn = float(1.0)
    with pytest.raises(TypeError):
        selcfg.xconfig_decn = str("s")
    with pytest.raises(TypeError):
        selcfg.xconfig_decn = numpy.random.random(common_nconfig)
    with pytest.raises(TypeError):
        selcfg.xconfig_decn = common_xconfig_decn_real

def test_xconfig_decn_fset_ValueError(selcfg, common_xconfig_decn_integer):
    with pytest.raises(ValueError):
        selcfg.xconfig_decn = common_xconfig_decn_integer

def test_xconfig_decn_fdel(selcfg):
    with pytest.raises(AttributeError):
        del selcfg.xconfig_decn

############################## Test class methods ##############################
def test_BinarySelectionConfiguration_sample_xconfig_is_concrete():
    assert_method_isconcrete(BinarySelectionConfiguration, "sample_xconfig")

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
