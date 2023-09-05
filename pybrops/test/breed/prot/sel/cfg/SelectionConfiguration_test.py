import pytest
from pybrops.test.assert_python import assert_concrete_class, assert_concrete_property, assert_docstring, assert_semiabstract_class, not_raises
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
def test_SelectionConfiguration_is_concrete():
    assert_concrete_class(SelectionConfiguration)

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

def test_ncross_fdel(selcfg):
    with pytest.raises(AttributeError):
        del selcfg.ncross

### nparent ###
def test_SelectionConfiguration_nparent_is_concrete():
    assert_concrete_property(SelectionConfiguration, "nparent")

def test_nparent_fget(selcfg, common_nparent):
    assert selcfg.nparent == common_nparent

def test_nparent_fset(selcfg, common_nparent):
    with not_raises(Exception):
        selcfg.nparent = common_nparent
    with not_raises(Exception):
        selcfg.nparent = int(1)
    with not_raises(Exception):
        selcfg.nparent = numpy.int8(1)
    with not_raises(Exception):
        selcfg.nparent = numpy.int16(1)
    with not_raises(Exception):
        selcfg.nparent = numpy.int32(1)
    with not_raises(Exception):
        selcfg.nparent = numpy.int64(1)

def test_nparent_fset_TypeError(selcfg):
    with pytest.raises(TypeError):
        selcfg.nparent = object()
    with pytest.raises(TypeError):
        selcfg.nparent = None
    with pytest.raises(TypeError):
        selcfg.nparent = float(1.0)
    with pytest.raises(TypeError):
        selcfg.nparent = str("s")

def test_nparent_fset_ValueError(selcfg):
    with pytest.raises(ValueError):
        selcfg.nparent = int(-1)
    with pytest.raises(ValueError):
        selcfg.nparent = int(0)

def test_nparent_fdel(selcfg):
    with pytest.raises(AttributeError):
        del selcfg.nparent

### nmating ###
def test_SelectionConfiguration_nmating_is_concrete():
    assert_concrete_property(SelectionConfiguration, "nmating")

def test_nmating_fget(selcfg, common_nmating):
    assert numpy.all(selcfg.nmating == common_nmating)

def test_nmating_fset(selcfg, common_nmating, common_ncross):
    with not_raises(Exception):
        selcfg.nmating = common_nmating
    with not_raises(Exception):
        selcfg.nmating = int(1)
    with not_raises(Exception):
        selcfg.nmating = numpy.int8(1)
    with not_raises(Exception):
        selcfg.nmating = numpy.int16(1)
    with not_raises(Exception):
        selcfg.nmating = numpy.int32(1)
    with not_raises(Exception):
        selcfg.nmating = numpy.int64(1)
    with not_raises(Exception):
        selcfg.nmating = numpy.repeat(numpy.int8(1), common_ncross)
    with not_raises(Exception):
        selcfg.nmating = numpy.repeat(numpy.int16(1), common_ncross)
    with not_raises(Exception):
        selcfg.nmating = numpy.repeat(numpy.int32(1), common_ncross)
    with not_raises(Exception):
        selcfg.nmating = numpy.repeat(numpy.int64(1), common_ncross)

def test_nmating_fset_TypeError(selcfg, common_ncross):
    with pytest.raises(TypeError):
        selcfg.nmating = object()
    with pytest.raises(TypeError):
        selcfg.nmating = None
    with pytest.raises(TypeError):
        selcfg.nmating = float(1.0)
    with pytest.raises(TypeError):
        selcfg.nmating = str("s")
    with pytest.raises(TypeError):
        selcfg.nmating = numpy.repeat(float(1.0), common_ncross)

def test_nmating_fset_ValueError(selcfg, common_nmating, common_ncross):
    with pytest.raises(ValueError):
        selcfg.nmating = int(-1)
    with pytest.raises(ValueError):
        selcfg.nmating = int(0)
    with pytest.raises(ValueError):
        selcfg.nmating = numpy.repeat(common_nmating, common_ncross+1)
    with pytest.raises(ValueError):
        selcfg.nmating = numpy.repeat(int(-1), common_ncross)
    with pytest.raises(ValueError):
        selcfg.nmating = numpy.repeat(int(0), common_ncross)

def test_nmating_fdel(selcfg):
    with pytest.raises(AttributeError):
        del selcfg.nmating

### nprogeny ###
def test_SelectionConfiguration_nprogeny_is_concrete():
    assert_concrete_property(SelectionConfiguration, "nprogeny")

def test_nprogeny_fget(selcfg, common_nprogeny):
    assert numpy.all(selcfg.nprogeny == common_nprogeny)

def test_nprogeny_fset(selcfg, common_nprogeny, common_ncross):
    with not_raises(Exception):
        selcfg.nprogeny = common_nprogeny
    with not_raises(Exception):
        selcfg.nprogeny = int(1)
    with not_raises(Exception):
        selcfg.nprogeny = numpy.int8(1)
    with not_raises(Exception):
        selcfg.nprogeny = numpy.int16(1)
    with not_raises(Exception):
        selcfg.nprogeny = numpy.int32(1)
    with not_raises(Exception):
        selcfg.nprogeny = numpy.int64(1)
    with not_raises(Exception):
        selcfg.nprogeny = numpy.repeat(numpy.int8(1), common_ncross)
    with not_raises(Exception):
        selcfg.nprogeny = numpy.repeat(numpy.int16(1), common_ncross)
    with not_raises(Exception):
        selcfg.nprogeny = numpy.repeat(numpy.int32(1), common_ncross)
    with not_raises(Exception):
        selcfg.nprogeny = numpy.repeat(numpy.int64(1), common_ncross)

def test_nprogeny_fset_TypeError(selcfg, common_ncross):
    with pytest.raises(TypeError):
        selcfg.nprogeny = object()
    with pytest.raises(TypeError):
        selcfg.nprogeny = None
    with pytest.raises(TypeError):
        selcfg.nprogeny = float(1.0)
    with pytest.raises(TypeError):
        selcfg.nprogeny = str("s")
    with pytest.raises(TypeError):
        selcfg.nprogeny = numpy.repeat(float(1.0), common_ncross)

def test_nprogeny_fset_ValueError(selcfg, common_nprogeny, common_ncross):
    with pytest.raises(ValueError):
        selcfg.nprogeny = int(-1)
    with pytest.raises(ValueError):
        selcfg.nprogeny = int(0)
    with pytest.raises(ValueError):
        selcfg.nprogeny = numpy.repeat(common_nprogeny, common_ncross+1)
    with pytest.raises(ValueError):
        selcfg.nprogeny = numpy.repeat(int(-1), common_ncross)
    with pytest.raises(ValueError):
        selcfg.nprogeny = numpy.repeat(int(0), common_ncross)

def test_nprogeny_fdel(selcfg):
    with pytest.raises(AttributeError):
        del selcfg.nprogeny

### pgmat ###
def test_SelectionConfiguration_pgmat_is_concrete():
    assert_concrete_property(SelectionConfiguration, "pgmat")

def test_pgmat_fget(selcfg, common_pgmat):
    assert numpy.all(selcfg.pgmat == common_pgmat)

def test_pgmat_fset(selcfg, common_pgmat):
    with not_raises(Exception):
        selcfg.pgmat = common_pgmat

def test_pgmat_fset_TypeError(selcfg):
    with pytest.raises(TypeError):
        selcfg.pgmat = object()
    with pytest.raises(TypeError):
        selcfg.pgmat = None
    with pytest.raises(TypeError):
        selcfg.pgmat = float(1.0)
    with pytest.raises(TypeError):
        selcfg.pgmat = str("s")

def test_pgmat_fdel(selcfg):
    with pytest.raises(AttributeError):
        del selcfg.pgmat

### xconfig ###
def test_SelectionConfiguration_xconfig_is_concrete():
    assert_concrete_property(SelectionConfiguration, "xconfig")

def test_xconfig_fget(selcfg, common_xconfig):
    assert numpy.all(selcfg.xconfig == common_xconfig)

def test_xconfig_fset(selcfg, common_xconfig, common_ncross, common_nparent, common_ntaxa):
    with not_raises(Exception):
        selcfg.xconfig = common_xconfig
    with not_raises(Exception):
        selcfg.xconfig = numpy.random.randint(0, common_ntaxa, (common_ncross, common_nparent))

def test_xconfig_fset_TypeError(selcfg, common_ncross, common_nparent):
    with pytest.raises(TypeError):
        selcfg.xconfig = object()
    with pytest.raises(TypeError):
        selcfg.xconfig = None
    with pytest.raises(TypeError):
        selcfg.xconfig = float(1.0)
    with pytest.raises(TypeError):
        selcfg.xconfig = str("s")
    with pytest.raises(TypeError):
        selcfg.xconfig = numpy.random.random((common_ncross, common_nparent))

def test_xconfig_fset_ValueError(selcfg, common_ncross, common_nparent, common_ntaxa):
    with pytest.raises(ValueError):
        selcfg.xconfig = numpy.random.randint(0, common_ntaxa, (common_ncross+1, common_nparent+1))

def test_xconfig_fdel(selcfg):
    with pytest.raises(AttributeError):
        del selcfg.xconfig

############################# Test class utilities #############################
def test_check_is_SelectionConfiguration(selcfg):
    check_is_SelectionConfiguration(selcfg, "selcfg")
