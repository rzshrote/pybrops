import numpy
import pytest
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function
from pybrops.test.breed.prot.mate.common_fixtures import *
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.breed.prot.mate.FourWayCross import FourWayCross
from pybrops.breed.prot.mate.FourWayCross import check_is_FourWayCross

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mprot(common_rng):
    yield FourWayCross(rng = common_rng)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(FourWayCross)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(FourWayCross, "__init__")

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

def test_check_is_FourWayCross_is_concrete():
    assert_concrete_function(check_is_FourWayCross)

def test_check_is_FourWayCross(mprot):
    with not_raises(TypeError):
        check_is_FourWayCross(mprot, "mprot")
    with pytest.raises(TypeError):
        check_is_FourWayCross(object(), "mprot")
    with pytest.raises(TypeError):
        check_is_FourWayCross(None, "mprot")
