import numpy
import pytest
from numpy.random import Generator, RandomState, PCG64
from pybrops.test.assert_python import assert_abstract_method, assert_concrete_property, assert_docstring, assert_mixin_class, not_raises
from pybrops.breed.prot.sel.cfg.SampledSelectionConfigurationMixin import SampledSelectionConfigurationMixin

from .common_fixtures import *

class DummySampledSelectionConfigurationMixin(SampledSelectionConfigurationMixin):
    def __init__(self, xconfig_decn, rng):
        self.xconfig_decn = xconfig_decn
        self.rng = rng
    def sample_xconfig(self, return_xconfig: bool):
        return super().sample_xconfig(return_xconfig)

@pytest.fixture
def selcfg(
        common_xconfig_decn_generic,
        common_rng
    ):
    out = DummySampledSelectionConfigurationMixin(
        xconfig_decn = common_xconfig_decn_generic,
        rng = common_rng
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_SampledSelectionConfigurationMixin_is_mixin():
    assert_mixin_class(SampledSelectionConfigurationMixin)

############################## Test class docstring ############################
def test_SampledSelectionConfigurationMixin_docstring():
    assert_docstring(SampledSelectionConfigurationMixin)

############################# Test class properties ############################

### xconfig_decn ###
def test_SampledSelectionConfigurationMixin_xconfig_decn_is_concrete():
    assert_concrete_property(SampledSelectionConfigurationMixin, "xconfig_decn")

def test_xconfig_decn_fget(selcfg, common_xconfig_decn_generic):
    assert numpy.all(selcfg.xconfig_decn == common_xconfig_decn_generic)

def test_xconfig_decn_fset(selcfg, common_xconfig_decn_generic, common_nconfig, common_ntaxa):
    with not_raises(Exception):
        selcfg.xconfig_decn = common_xconfig_decn_generic
    with not_raises(Exception):
        selcfg.xconfig_decn = numpy.random.randint(0, common_ntaxa, common_nconfig)

def test_xconfig_decn_fset_TypeError(selcfg):
    with pytest.raises(TypeError):
        selcfg.xconfig_decn = object()
    with pytest.raises(TypeError):
        selcfg.xconfig_decn = None
    with pytest.raises(TypeError):
        selcfg.xconfig_decn = float(1.0)
    with pytest.raises(TypeError):
        selcfg.xconfig_decn = str("s")

def test_xconfig_decn_fdel(selcfg):
    with pytest.raises(AttributeError):
        del selcfg.xconfig_decn

### rng ###
def test_SampledSelectionConfigurationMixin_rng_is_concrete():
    assert_concrete_property(SampledSelectionConfigurationMixin, "rng")

def test_rng_fget(selcfg):
    assert isinstance(selcfg.rng, Generator) or isinstance(selcfg.rng, RandomState)

def test_rng_fset(selcfg, common_rng):
    with not_raises(Exception):
        selcfg.rng = common_rng
    with not_raises(Exception):
        selcfg.rng = Generator(PCG64(1))
    with not_raises(Exception):
        selcfg.rng = RandomState(1)
    with not_raises(Exception):
        selcfg.rng = None

def test_rng_fset_TypeError(selcfg):
    with pytest.raises(TypeError):
        selcfg.rng = object()
    with pytest.raises(TypeError):
        selcfg.rng = float(1.0)
    with pytest.raises(TypeError):
        selcfg.rng = str("s")

def test_rng_fdel(selcfg):
    with pytest.raises(AttributeError):
        del selcfg.rng

############################## Test class methods ##############################
def test_SampledSelectionConfigurationMixin_sample_xconfig_is_abstract():
    assert_abstract_method(SampledSelectionConfigurationMixin, "sample_xconfig")
