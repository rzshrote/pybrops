import numpy
import pytest
from numpy.random import PCG64
from numpy.random import Generator

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.breed.prot.mate.FourWayDHCross import FourWayDHCross
from pybrops.breed.prot.mate.FourWayDHCross import is_FourWayDHCross
from pybrops.breed.prot.mate.FourWayDHCross import check_is_FourWayDHCross

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
def test_mate(mprot, dpgvmat, sel, rng):
    progeny = mprot.mate(dpgvmat, sel, 1, 2, nself = 0)
    # print("parents:\n", dpgvmat.mat)
    # print("progeny:\n", progeny.mat)
    # raise RuntimeError("stop")
    assert isinstance(progeny, DensePhasedGenotypeMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_FourWayDHCross_is_concrete():
    assert_concrete_function(is_FourWayDHCross)

def test_is_FourWayDHCross(mprot):
    assert is_FourWayDHCross(mprot)

def test_check_is_FourWayDHCross_is_concrete():
    assert_concrete_function(check_is_FourWayDHCross)

def test_check_is_FourWayDHCross(mprot):
    with not_raises(TypeError):
        check_is_FourWayDHCross(mprot, "mprot")
    with pytest.raises(TypeError):
        check_is_FourWayDHCross(None, "mprot")
