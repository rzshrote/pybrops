import pytest
from pybrops.test.assert_python import assert_concrete_class, assert_concrete_property, assert_docstring, assert_semiabstract_class, not_raises
from pybrops.breed.prot.sel.cfg.MateSelectionConfiguration import MateSelectionConfiguration
from pybrops.breed.prot.sel.cfg.MateSelectionConfiguration import check_is_MateSelectionConfiguration

from pybrops.test.breed.prot.sel.cfg.common_fixtures import *

class DummyMateSelectionConfiguration(MateSelectionConfiguration):
    def __init__(self, ncross, nparent, nmating, nprogeny, pgmat, xconfig, xconfig_xmap):
        self.ncross = ncross
        self.nparent = nparent
        self.nmating = nmating
        self.nprogeny = nprogeny
        self.pgmat = pgmat
        self.xconfig = xconfig
        self.xconfig_xmap = xconfig_xmap

@pytest.fixture
def selcfg(
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_pgmat,
        common_xconfig,
        common_xconfig_xmap
    ):
    out = DummyMateSelectionConfiguration(
        ncross = common_ncross,
        nparent = common_nparent,
        nmating = common_nmating,
        nprogeny = common_nprogeny,
        pgmat = common_pgmat,
        xconfig = common_xconfig,
        xconfig_xmap = common_xconfig_xmap
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_MateSelectionConfiguration_is_concrete():
    assert_concrete_class(MateSelectionConfiguration)

############################## Test class docstring ############################
def test_MateSelectionConfiguration_docstring():
    assert_docstring(MateSelectionConfiguration)

############################# Test class properties ############################

### xconfig_xmap ###
def test_MateSelectionConfiguration_xconfig_xmap_is_concrete():
    assert_concrete_property(MateSelectionConfiguration, "xconfig_xmap")

def test_xconfig_xmap_fget(selcfg, common_xconfig_xmap):
    assert numpy.all(selcfg.xconfig_xmap == common_xconfig_xmap)

def test_xconfig_xmap_fset(selcfg, common_xconfig_xmap, common_nconfig, common_nparent, common_ntaxa):
    with not_raises(Exception):
        selcfg.xconfig_xmap = common_xconfig_xmap
    with not_raises(Exception):
        selcfg.xconfig_xmap = numpy.random.randint(0, common_ntaxa, (common_nconfig, common_nparent))

def test_xconfig_xmap_fset_TypeError(selcfg, common_nconfig, common_nparent):
    with pytest.raises(TypeError):
        selcfg.xconfig_xmap = object()
    with pytest.raises(TypeError):
        selcfg.xconfig_xmap = None
    with pytest.raises(TypeError):
        selcfg.xconfig_xmap = float(1.0)
    with pytest.raises(TypeError):
        selcfg.xconfig_xmap = str("s")
    with pytest.raises(TypeError):
        selcfg.xconfig_xmap = numpy.random.random((common_nconfig, common_nparent))

def test_xconfig_xmap_fset_ValueError(selcfg, common_nconfig, common_nparent, common_ntaxa):
    with pytest.raises(ValueError):
        selcfg.xconfig_xmap = numpy.random.randint(0, common_ntaxa, (common_nconfig+1, common_nparent+1))

def test_xconfig_xmap_fdel(selcfg):
    with pytest.raises(AttributeError):
        del selcfg.xconfig_xmap

############################# Test class utilities #############################
def test_check_is_MateSelectionConfiguration(selcfg):
    check_is_MateSelectionConfiguration(selcfg, "selcfg")
