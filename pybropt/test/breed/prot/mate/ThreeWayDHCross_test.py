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
from pybropt.popgen.gmap import ExtendedGeneticMap
from pybropt.popgen.gmap import HaldaneMapFunction
from pybropt.breed.prot.mate import ThreeWayDHCross
from pybropt.breed.prot.mate import is_ThreeWayDHCross
from pybropt.breed.prot.mate import check_is_ThreeWayDHCross
from pybropt.breed.prot.mate import cond_check_is_ThreeWayDHCross

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
    yield ThreeWayDHCross(
        rng = rng
    )

@pytest.fixture
def sel():
    yield numpy.int64([0,1,2,0,2,2,1,2,1])

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(ThreeWayDHCross)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(ThreeWayDHCross, "__init__")

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
    progeny, misc = mprot.mate(dpgvmat, sel, 1, 2, s = 0)
    # print("parents:\n", dpgvmat.mat)
    # print("progeny:\n", progeny.mat)
    # raise RuntimeError("stop")
    assert is_DensePhasedGenotypeMatrix(progeny)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_ThreeWayDHCross_is_concrete():
    generic_assert_concrete_function(is_ThreeWayDHCross)

def test_is_ThreeWayDHCross(mprot):
    assert is_ThreeWayDHCross(mprot)

def test_check_is_ThreeWayDHCross_is_concrete():
    generic_assert_concrete_function(check_is_ThreeWayDHCross)

def test_check_is_ThreeWayDHCross(mprot):
    with not_raises(TypeError):
        check_is_ThreeWayDHCross(mprot, "mprot")
    with pytest.raises(TypeError):
        check_is_ThreeWayDHCross(None, "mprot")

def test_cond_check_is_ThreeWayDHCross_is_concrete():
    generic_assert_concrete_function(cond_check_is_ThreeWayDHCross)
