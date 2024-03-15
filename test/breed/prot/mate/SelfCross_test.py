import numpy
import pytest
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete
from .common_fixtures import *
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.breed.prot.mate.SelfCross import SelfCross
from pybrops.breed.prot.mate.SelfCross import check_is_SelfCross

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mprot(common_rng):
    yield SelfCross(rng = common_rng)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(SelfCross)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test___init___is_concrete():
    assert_method_isconcrete(SelfCross, "__init__")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_mate(mprot, common_dpgmat, common_xconfig_self):
    # parameters
    nmating = 2
    nprogeny = 2
    nself = 2

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_self, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )

    assert isinstance(progeny, DensePhasedGenotypeMatrix)

def test_mate_nmating(mprot, common_dpgmat, common_xconfig_self):
    # parameters
    nmating = 10
    nprogeny = 1
    nself = 0

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_self, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )
    assert progeny.ntaxa == len(common_xconfig_self) * nmating * nprogeny

def test_mate_nmating_nself(mprot, common_dpgmat, common_xconfig_self):
    # parameters
    nmating = 10
    nprogeny = 1
    nself = 1

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_self, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )
    assert progeny.ntaxa == len(common_xconfig_self) * nmating * nprogeny

def test_mate_nprogeny(mprot, common_dpgmat, common_xconfig_self):
    # parameters
    nmating = 1
    nprogeny = 10
    nself = 0

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_self, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )
    assert progeny.ntaxa == len(common_xconfig_self) * nmating * nprogeny

def test_mate_nprogeny_nself(mprot, common_dpgmat, common_xconfig_self):
    # parameters
    nmating = 1
    nprogeny = 10
    nself = 1

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_self, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )

    assert progeny.ntaxa == len(common_xconfig_self) * nmating * nprogeny

def test_mate_nmating_nprogeny(mprot, common_dpgmat, common_xconfig_self):
    # parameters
    nmating = 10
    nprogeny = 10
    nself = 0

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_self, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )

    assert progeny.ntaxa == len(common_xconfig_self) * nmating * nprogeny

def test_mate_nmating_nprogeny_nself(mprot, common_dpgmat, common_xconfig_self):
    # parameters
    nmating = 10
    nprogeny = 10
    nself = 1

    # make progenies
    progeny = mprot.mate(
        pgmat = common_dpgmat, 
        xconfig = common_xconfig_self, 
        nmating = nmating, 
        nprogeny = nprogeny, 
        miscout = None, 
        nself = nself
    )

    assert progeny.ntaxa == len(common_xconfig_self) * nmating * nprogeny

################################################################################
######################### Test class utility functions #########################
################################################################################

def test_check_is_SelfCross_is_concrete():
    assert_function_isconcrete(check_is_SelfCross)

def test_check_is_SelfCross(mprot):
    with not_raises(TypeError):
        check_is_SelfCross(mprot, "mprot")
    with pytest.raises(TypeError):
        check_is_SelfCross(object(), "mprot")
    with pytest.raises(TypeError):
        check_is_SelfCross(None, "mprot")
