import numpy
import pytest
from numpy.random import PCG64
from numpy.random import Generator

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.popgen.gmat import DensePhasedGenotypeMatrix
from pybropt.popgen.gmat import is_DensePhasedGenotypeMatrix
from pybropt.popgen.gmap import ExtendedGeneticMap
from pybropt.popgen.gmap import HaldaneMapFunction
from pybropt.breed.prot.mate import TwoWayDHCross
from pybropt.breed.prot.mate import is_TwoWayDHCross
from pybropt.breed.prot.mate import check_is_TwoWayDHCross
from pybropt.breed.prot.mate import cond_check_is_TwoWayDHCross

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def egmap(shared_datadir):
    e = ExtendedGeneticMap.from_egmap(shared_datadir / "sample.egmap")
    e.group()
    e.build_spline(kind = 'linear', fill_value = 'extrapolate')
    yield e

@pytest.fixture
def gmapfn():
    yield HaldaneMapFunction()

@pytest.fixture
def dpgvmat(shared_datadir, egmap, gmapfn):
    data_path = shared_datadir / "sample.vcf"
    dpgvmat = DensePhasedGenotypeMatrix.from_vcf(data_path)
    dpgvmat.group()
    dpgvmat.interp_xoprob(egmap, gmapfn)
    yield dpgvmat

@pytest.fixture
def rng():
    yield Generator(PCG64(543212345))

@pytest.fixture
def mprot(rng):
    yield TwoWayDHCross(
        rng = rng
    )

@pytest.fixture
def sel():
    yield numpy.int64([0,1,0,2,1,2])

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(TwoWayDHCross)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(TwoWayDHCross, "__init__")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_mate(mprot, dpgvmat, sel, rng):
    ncross = 2
    nprogeny = 2
    s = 2
    progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, s)
    # print("parents:\n", dpgvmat.mat)
    # print("progeny:\n", progeny.mat)
    # raise RuntimeError("stop")
    assert is_DensePhasedGenotypeMatrix(progeny)
    mat = progeny.mat
    for i in range(1,len(mat)):
        assert numpy.all(mat[0] == mat[i])

def test_mate_ncross(mprot, dpgvmat, sel, rng):
    ncross = 10
    nprogeny = 1
    s = 0
    progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, s)
    assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

def test_mate_ncross_s(mprot, dpgvmat, sel, rng):
    ncross = 10
    nprogeny = 1
    s = 1
    progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, s)
    assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

def test_mate_nprogeny(mprot, dpgvmat, sel, rng):
    ncross = 1
    nprogeny = 10
    s = 0
    progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, s)
    assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

def test_mate_nprogeny_s(mprot, dpgvmat, sel, rng):
    ncross = 1
    nprogeny = 10
    s = 1
    progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, s)
    assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

def test_mate_ncross_nprogeny(mprot, dpgvmat, sel, rng):
    ncross = 10
    nprogeny = 10
    s = 0
    progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, s)
    assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

def test_mate_ncross_nprogeny_s(mprot, dpgvmat, sel, rng):
    ncross = 10
    nprogeny = 10
    s = 1
    progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, s)
    assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_TwoWayDHCross_is_concrete():
    generic_assert_concrete_function(is_TwoWayDHCross)

def test_is_TwoWayDHCross(mprot):
    assert is_TwoWayDHCross(mprot)

def test_check_is_TwoWayDHCross_is_concrete():
    generic_assert_concrete_function(check_is_TwoWayDHCross)

def test_check_is_TwoWayDHCross(mprot):
    with not_raises(TypeError):
        check_is_TwoWayDHCross(mprot, "mprot")
    with pytest.raises(TypeError):
        check_is_TwoWayDHCross(None, "mprot")

def test_cond_check_is_TwoWayDHCross_is_concrete():
    generic_assert_concrete_function(cond_check_is_TwoWayDHCross)
