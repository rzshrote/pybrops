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
from pybrops.breed.prot.mate.TwoWayCross import TwoWayCross
from pybrops.breed.prot.mate.TwoWayCross import check_is_TwoWayCross

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
    yield TwoWayCross(
        rng = rng
    )

@pytest.fixture
def sel():
    yield numpy.int64([0,1,0,2,1,2])

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(TwoWayCross)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(TwoWayCross, "__init__")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_mate(mprot, dpgvmat, sel):
    ncross = 2
    nprogeny = 2
    nself = 2
    progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, nself)
    # print("parents:\n", dpgvmat.mat)
    # print("progeny:\n", progeny.mat)
    # raise RuntimeError("stop")
    assert isinstance(progeny, DensePhasedGenotypeMatrix)
    mat = progeny.mat
    for i in range(1,len(mat)):
        assert numpy.any(mat[0] != mat[i])

# def test_mate_ncross(mprot, dpgvmat, sel, rng):
#     ncross = 10
#     nprogeny = 1
#     nself = 0
#     progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, nself)
#     assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

# def test_mate_ncross_s(mprot, dpgvmat, sel, rng):
#     ncross = 10
#     nprogeny = 1
#     nself = 1
#     progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, nself)
#     assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

# def test_mate_nprogeny(mprot, dpgvmat, sel, rng):
#     ncross = 1
#     nprogeny = 10
#     nself = 0
#     progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, nself)
#     assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

# def test_mate_nprogeny_s(mprot, dpgvmat, sel, rng):
#     ncross = 1
#     nprogeny = 10
#     nself = 1
#     progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, nself)
#     assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

# def test_mate_ncross_nprogeny(mprot, dpgvmat, sel, rng):
#     ncross = 10
#     nprogeny = 10
#     nself = 0
#     progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, nself)
#     assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

# def test_mate_ncross_nprogeny_s(mprot, dpgvmat, sel, rng):
#     ncross = 10
#     nprogeny = 10
#     nself = 1
#     progeny = mprot.mate(dpgvmat, sel, ncross, nprogeny, nself)
#     assert progeny.ntaxa == (len(sel) // 2) * ncross * nprogeny

# ################################################################################
# ######################### Test class utility functions #########################
# ################################################################################

# def test_check_is_TwoWayCross_is_concrete():
#     generic_assert_concrete_function(check_is_TwoWayCross)

# def test_check_is_TwoWayCross(mprot):
#     with not_raises(TypeError):
#         check_is_TwoWayCross(mprot, "mprot")
#     with pytest.raises(TypeError):
#         check_is_TwoWayCross(None, "mprot")
