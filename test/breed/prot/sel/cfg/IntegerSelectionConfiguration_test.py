import numpy
import pytest
from pybrops.test.assert_python import assert_method_isconcrete, assert_module_documentation, assert_property_isconcrete, assert_class_documentation, assert_class_isconcrete, not_raises
from pybrops.breed.prot.sel.cfg.IntegerSelectionConfiguration import IntegerSelectionConfiguration

from .common_fixtures import *

@pytest.fixture
def selcfg(
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_pgmat,
        common_xconfig_decn_integer,
        common_rng
    ):
    out = IntegerSelectionConfiguration(
        ncross = common_ncross,
        nparent = common_nparent,
        nmating = common_nmating,
        nprogeny = common_nprogeny,
        pgmat = common_pgmat,
        xconfig_decn = common_xconfig_decn_integer,
        rng = common_rng
    )
    yield out

############################ Test module attributes ############################
def test_IntegerSelectionConfiguration_module_is_concrete():
    import pybrops.breed.prot.sel.cfg.IntegerSelectionConfiguration
    assert_module_documentation(pybrops.breed.prot.sel.cfg.IntegerSelectionConfiguration)

############################ Test class attributes #############################
def test_IntegerSelectionConfiguration_is_concrete():
    assert_class_isconcrete(IntegerSelectionConfiguration)

##################### Test class special concrete methods ######################
def test_IntegerSelectionConfiguration___init___is_concrete():
    assert_method_isconcrete(IntegerSelectionConfiguration, "__init__")

############################# Test class properties ############################

### xconfig_decn ###
def test_IntegerSelectionConfiguration_xconfig_decn_is_concrete():
    assert_property_isconcrete(IntegerSelectionConfiguration, "xconfig_decn")

def test_xconfig_decn_fget(selcfg, common_xconfig_decn_integer):
    assert numpy.all(selcfg.xconfig_decn == common_xconfig_decn_integer)

def test_xconfig_decn_fset(selcfg, common_xconfig_decn_integer, common_xconfig_decn_binary, common_nconfig, common_ntaxa):
    with not_raises(Exception):
        selcfg.xconfig_decn = common_xconfig_decn_integer
    with not_raises(Exception):
        selcfg.xconfig_decn = common_xconfig_decn_binary
    with not_raises(Exception):
        selcfg.xconfig_decn = numpy.random.randint(0, common_ntaxa, common_nconfig)

def test_xconfig_decn_fset_TypeError(selcfg, common_xconfig_decn_bool, common_xconfig_decn_real, common_nconfig):
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
        selcfg.xconfig_decn = common_xconfig_decn_bool
    with pytest.raises(TypeError):
        selcfg.xconfig_decn = common_xconfig_decn_real

def test_xconfig_decn_fdel(selcfg):
    with pytest.raises(AttributeError):
        del selcfg.xconfig_decn

############################## Test class methods ##############################
def test_IntegerSelectionConfiguration_sample_xconfig_is_concrete():
    assert_method_isconcrete(IntegerSelectionConfiguration, "sample_xconfig")

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
