import numpy
import pytest
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function
from .common_fixtures import *
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.breed.prot.mate.FourWayDHCross import FourWayDHCross
from pybrops.breed.prot.mate.FourWayDHCross import check_is_FourWayDHCross

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mprot(common_rng):
    yield FourWayDHCross(rng = common_rng)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(FourWayDHCross)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(FourWayDHCross, "__init__")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_mate(mprot, common_dpgmat, common_xconfig_fourway):
    # parameters
    nmating = 2
    nprogeny = 2
    nself = 2

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_fourway, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )

    assert isinstance(progeny, DensePhasedGenotypeMatrix)

    mat = progeny.mat
    for i in range(1,len(mat)):
        assert numpy.all(mat[0] == mat[i])

def test_mate_nmating(mprot, common_dpgmat, common_xconfig_fourway):
    # parameters
    nmating = 10
    nprogeny = 1
    nself = 0

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_fourway, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )
    assert progeny.ntaxa == len(common_xconfig_fourway) * nmating * nprogeny

def test_mate_nmating_nself(mprot, common_dpgmat, common_xconfig_fourway):
    # parameters
    nmating = 10
    nprogeny = 1
    nself = 1

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_fourway, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )
    assert progeny.ntaxa == len(common_xconfig_fourway) * nmating * nprogeny

def test_mate_nprogeny(mprot, common_dpgmat, common_xconfig_fourway):
    # parameters
    nmating = 1
    nprogeny = 10
    nself = 0

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_fourway, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )
    assert progeny.ntaxa == len(common_xconfig_fourway) * nmating * nprogeny

def test_mate_nprogeny_nself(mprot, common_dpgmat, common_xconfig_fourway):
    # parameters
    nmating = 1
    nprogeny = 10
    nself = 1

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_fourway, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )

    assert progeny.ntaxa == len(common_xconfig_fourway) * nmating * nprogeny

def test_mate_nmating_nprogeny(mprot, common_dpgmat, common_xconfig_fourway):
    # parameters
    nmating = 10
    nprogeny = 10
    nself = 0

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_fourway, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )

    assert progeny.ntaxa == len(common_xconfig_fourway) * nmating * nprogeny

def test_mate_nmating_nprogeny_nself(mprot, common_dpgmat, common_xconfig_fourway):
    # parameters
    nmating = 10
    nprogeny = 10
    nself = 1

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_fourway, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )

    assert progeny.ntaxa == len(common_xconfig_fourway) * nmating * nprogeny

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_FourWayDHCross_is_concrete():
    assert_concrete_function(check_is_FourWayDHCross)

def test_check_is_FourWayDHCross(mprot):
    with not_raises(TypeError):
        check_is_FourWayDHCross(mprot, "mprot")
    with pytest.raises(TypeError):
        check_is_FourWayDHCross(object(), "mprot")
    with pytest.raises(TypeError):
        check_is_FourWayDHCross(None, "mprot")