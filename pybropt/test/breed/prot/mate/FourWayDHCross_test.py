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
from pybropt.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap
from pybropt.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybropt.breed.prot.mate.FourWayDHCross import FourWayDHCross
from pybropt.breed.prot.mate.FourWayDHCross import is_FourWayDHCross
from pybropt.breed.prot.mate.FourWayDHCross import check_is_FourWayDHCross
from pybropt.breed.prot.mate.FourWayDHCross import cond_check_is_FourWayDHCross

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
    yield FourWayDHCross(
        rng = rng
    )

@pytest.fixture
def sel():
    yield numpy.int64([0,1,2,1, 0,2,2,1, 1,2,1,0])

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(FourWayDHCross)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(FourWayDHCross, "__init__")

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
    progeny = mprot.mate(dpgvmat, sel, 1, 2, s = 0)
    # print("parents:\n", dpgvmat.mat)
    # print("progeny:\n", progeny.mat)
    # raise RuntimeError("stop")
    assert is_DensePhasedGenotypeMatrix(progeny)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_FourWayDHCross_is_concrete():
    generic_assert_concrete_function(is_FourWayDHCross)

def test_is_FourWayDHCross(mprot):
    assert is_FourWayDHCross(mprot)

def test_check_is_FourWayDHCross_is_concrete():
    generic_assert_concrete_function(check_is_FourWayDHCross)

def test_check_is_FourWayDHCross(mprot):
    with not_raises(TypeError):
        check_is_FourWayDHCross(mprot, "mprot")
    with pytest.raises(TypeError):
        check_is_FourWayDHCross(None, "mprot")

def test_cond_check_is_FourWayDHCross_is_concrete():
    generic_assert_concrete_function(cond_check_is_FourWayDHCross)
