import pytest
from pybrops.test.assert_python import assert_concrete_property, assert_docstring, assert_semiabstract_class, not_raises
from pybrops.breed.prot.sel.cfg.SelectionConfiguration import SelectionConfiguration
from pybrops.breed.prot.sel.cfg.SelectionConfiguration import check_is_SelectionConfiguration

from pybrops.test.breed.prot.sel.cfg.common_fixtures import *

class DummySelectionConfiguration(SelectionConfiguration):
    def __init__(self, ncross, nparent, nmating, nprogeny, pgmat, xconfig):
        self.ncross = ncross
        self.nparent = nparent
        self.nmating = nmating
        self.nprogeny = nprogeny
        self.pgmat = pgmat
        self.xconfig = xconfig

@pytest.fixture
def selcfg(
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_pgmat,
        common_xconfig
    ):
    out = DummySelectionConfiguration(
        ncross = common_ncross,
        nparent = common_nparent,
        nmating = common_nmating,
        nprogeny = common_nprogeny,
        pgmat = common_pgmat,
        xconfig = common_xconfig
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_SelectionConfiguration_is_semiabstract():
    assert_semiabstract_class(SelectionConfiguration)

############################## Test class docstring ############################
def test_SelectionConfiguration_docstring():
    assert_docstring(SelectionConfiguration)

############################# Test class properties ############################

### ncross ###
def test_SelectionConfiguration_ncross_is_concrete():
    assert_concrete_property(SelectionConfiguration, "ncross")

def test_ncross_fget(selcfg, common_ncross):
    assert selcfg.ncross == common_ncross

def test_ncross_fset(selcfg, common_ncross):
    with not_raises(Exception):
        selcfg.ncross = common_ncross
    with not_raises(Exception):
        selcfg.ncross = int(1)
    with not_raises(Exception):
        selcfg.ncross = numpy.int8(1)
    with not_raises(Exception):
        selcfg.ncross = numpy.int16(1)
    with not_raises(Exception):
        selcfg.ncross = numpy.int32(1)
    with not_raises(Exception):
        selcfg.ncross = numpy.int64(1)

def test_ncross_fset_TypeError(selcfg):
    with pytest.raises(TypeError):
        selcfg.ncross = object()
    with pytest.raises(TypeError):
        selcfg.ncross = None
    with pytest.raises(TypeError):
        selcfg.ncross = float(1.0)
    with pytest.raises(TypeError):
        selcfg.ncross = str("s")

def test_ncross_fset_ValueError(selcfg):
    with pytest.raises(ValueError):
        selcfg.ncross = int(-1)
    with pytest.raises(ValueError):
        selcfg.ncross = int(0)

### nparent ###
def test_SelectionConfiguration_nparent_is_concrete():
    assert_concrete_property(SelectionConfiguration, "nparent")

### nmating ###
def test_SelectionConfiguration_nmating_is_concrete():
    assert_concrete_property(SelectionConfiguration, "nmating")

### nprogeny ###
def test_SelectionConfiguration_nprogeny_is_concrete():
    assert_concrete_property(SelectionConfiguration, "nprogeny")


############################# Test class utilities #############################
def test_check_is_SelectionConfiguration(selcfg):
    check_is_SelectionConfiguration(selcfg, "selcfg")
