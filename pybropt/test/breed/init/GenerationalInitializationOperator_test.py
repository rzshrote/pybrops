import numpy
import pytest

from numpy.random import Generator
from numpy.random import PCG64

from pybropt.popgen.gmap import ExtendedGeneticMap
from pybropt.popgen.gmap import HaldaneMapFunction
from pybropt.popgen.gmat import DensePhasedGenotypeVariantMatrix
from pybropt.model.gmod import GenericLinearGenomicModel
from pybropt.breed.sel import ConventionalPhenotypicParentSelection
from pybropt.breed.mate import GenerationalTwoWayDHCross
from pybropt.breed.intg import GenerationalGenotypeIntegrationOperator
from pybropt.breed.eval import NoGxEEvaluationOperator
from pybropt.breed.intg import GenerationalBreedingValueIntegrationOperator
from pybropt.breed.calibr import TrueGenomicModelCalibrationOperator
from pybropt.breed.sel import FamilyPhenotypicSurvivorSelection
from pybropt.breed.init import GenerationalInitializationOperator

################################################################################
################# GenerationalInitializationOperator fixtures ##################
################################################################################
@pytest.fixture
def gmap(shared_datadir):
    data_path = shared_datadir / "Song_2016.linear.M.egmap"
    result = ExtendedGeneticMap.from_egmap(data_path)
    result.group()
    result.build_spline()
    yield result

@pytest.fixture
def gmapfn():
    yield HaldaneMapFunction()

@pytest.fixture
def dpgvmat(shared_datadir, gmap, gmapfn):
    data_path = shared_datadir / "Song_2016_phased_chr_1000.vcf"
    result = DensePhasedGenotypeVariantMatrix.from_vcf(data_path)
    result.group()
    result.interp_xoprob(gmap, gmapfn)
    yield result

@pytest.fixture
def k_p():
    yield 40

@pytest.fixture
def ncross():
    yield 1

@pytest.fixture
def nprogeny():
    yield 80

@pytest.fixture
def gqlen():
    yield 3

@pytest.fixture
def gwind():
    yield 3

@pytest.fixture
def gmult():
    yield 1000

@pytest.fixture
def rng():
    yield Generator(PCG64(192837465))

@pytest.fixture
def gmod_true(rng):
    mu = rng.uniform(100, 200, (3,1))
    beta = rng.normal(0.0, 1.0, (1000,3))
    trait = numpy.object_(["yield", "protein", "oil"])
    model_name = "test_true"
    params = {}
    yield GenericLinearGenomicModel(
        mu = mu,
        beta = beta,
        trait = trait,
        model_name = model_name,
        params = params
    )

@pytest.fixture
def burnin():
    yield 20

@pytest.fixture
def t_max():
    yield 20

@pytest.fixture
def pselop(k_p, ncross, nprogeny, rng):
    yield ConventionalPhenotypicParentSelection(
        k_p = k_p,
        traitwt_p = numpy.float64([1.0, 1.0, 1.0]),
        ncross = ncross,
        nprogeny = nprogeny,
        rng = rng
    )

@pytest.fixture
def mateop(rng, gmult):
    yield GenerationalTwoWayDHCross(
        rng = rng,
        gmult = gmult
    )

@pytest.fixture
def gintgop(gqlen, gwind, gmult):
    yield GenerationalGenotypeIntegrationOperator(
        gqlen = gqlen,
        gwind = gwind,
        gmult = gmult
    )

@pytest.fixture
def evalop(dpgvmat, gmod_true, rng):
    yield NoGxEEvaluationOperator.from_h2(
        gmat = dpgvmat,
        lgmod = gmod_true,
        nenv = 4,
        h2 = 0.4,
        rng = rng
    )

@pytest.fixture
def bvintgop(gqlen, gwind, gmult):
    yield GenerationalBreedingValueIntegrationOperator(
        gqlen = gqlen,
        gwind = gwind,
        gmult = gmult
    )

@pytest.fixture
def calop():
    yield TrueGenomicModelCalibrationOperator()

@pytest.fixture
def sselop(rng):
    yield FamilyPhenotypicSurvivorSelection(
        k_f = 4,
        traitwt_f = numpy.float64([1.0, 1.0, 1.0]),
        rng = rng,
    )

@pytest.fixture
def seed_nsel(k_p):
    yield k_p

@pytest.fixture
def seed_ncross(ncross):
    yield ncross

@pytest.fixture
def seed_nprogeny(nprogeny):
    yield nprogeny

@pytest.fixture
def seed_ncross():
    yield 1

@pytest.fixture
def initop(dpgvmat, rng, seed_nsel, seed_ncross, seed_nprogeny, gqlen, gwind, gmod_true, burnin, t_max, pselop, mateop, gintgop, evalop, bvintgop, calop, sselop):
    result = GenerationalInitializationOperator.from_dpgvmat(
        dpgvmat = dpgvmat,
        rng = rng,
        seed_nsel = seed_nsel,
        seed_ncross = seed_ncross,
        seed_nprogeny = seed_nprogeny,
        gqlen = gqlen,
        gwind = gwind,
        gmod_true = gmod_true,
        burnin = burnin,
        t_max = t_max,
        pselop = pselop,
        mateop = mateop,
        gintgop = gintgop,
        evalop = evalop,
        bvintgop = bvintgop,
        calop = calop,
        sselop = sselop,
        replace = False
    )
    yield result

################################################################################
#################################### Tests #####################################
################################################################################
def test_initialize(initop):
    # potential problem with diversity for this test set
    initop.initialize()
    raise RuntimeError
