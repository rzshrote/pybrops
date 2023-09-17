import pytest
import numpy
from pybrops.core.random.prng import global_prng
from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

@pytest.fixture
def common_egmap(shared_datadir):
    e = ExtendedGeneticMap.from_egmap(shared_datadir / "sample.egmap")
    e.group()
    e.build_spline(kind = 'linear', fill_value = 'extrapolate')
    yield e

@pytest.fixture
def common_gmapfn():
    yield HaldaneMapFunction()

@pytest.fixture
def common_dpgmat(shared_datadir, common_egmap, common_gmapfn):
    data_path = shared_datadir / "sample.vcf"
    dpgmat = DensePhasedGenotypeMatrix.from_vcf(data_path)
    dpgmat.group()
    dpgmat.interp_xoprob(common_egmap, common_gmapfn)
    yield dpgmat

@pytest.fixture
def common_rng():
    yield global_prng

@pytest.fixture
def common_xconfig_self():
    yield numpy.int64([
        [0],
        [1],
        [0],
        [2],
        [1],
        [2]
    ])

@pytest.fixture
def common_xconfig_twoway():
    yield numpy.int64([
        [0,1],
        [0,2],
        [1,2]
    ])

@pytest.fixture
def common_xconfig_threeway():
    yield numpy.int64([
        [0,1,2],
        [0,2,2],
        [1,2,1]
    ])

@pytest.fixture
def common_xconfig_fourway():
    yield numpy.int64([
        [0,1,2,1],
        [0,2,2,1],
        [1,2,1,0]
    ])
