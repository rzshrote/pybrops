import numpy
import pytest

from numpy.random import PCG64
from numpy.random import Generator

from pybropt.popgen.gmat import DensePhasedGenotypeMatrix
from pybropt.popgen.gmat import is_DensePhasedGenotypeMatrix

from pybropt.popgen.gmap import ExtendedGeneticMap
from pybropt.popgen.gmap import HaldaneMapFunction

from pybropt.breed.mate import GenerationalTwoWayDHCross

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
def twoway(rng):
    yield GenerationalTwoWayDHCross(
        rng = rng,
        gmult = 100
    )

@pytest.fixture
def sel():
    yield numpy.int64([0,1,0,2,1,2])

def test_mate(twoway, dpgvmat, sel, rng):
    progeny, misc = twoway.mate(1, 10, dpgvmat, sel, 1, 2, s = 0)
    # print(progeny.taxa_grp)
    # print("parents:\n", dpgvmat.mat)
    # print("progeny:\n", progeny.mat)
    # raise RuntimeError("stop")
