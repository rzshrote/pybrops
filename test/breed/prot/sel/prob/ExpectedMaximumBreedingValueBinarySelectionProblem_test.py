import numpy
import pytest

from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.breed.prot.mate.TwoWayDHCross import TwoWayDHCross
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_property_isconcrete
from pybrops.breed.prot.sel.prob.ExpectedMaximumBreedingValueSelectionProblem import ExpectedMaximumBreedingValueBinarySelectionProblem

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def ntaxa():
    yield 10

@pytest.fixture
def ntrait():
    yield 2

@pytest.fixture
def nvrnt():
    yield 1000

@pytest.fixture
def nparent():
    yield 2

@pytest.fixture
def ncross():
    yield 1

@pytest.fixture
def nprogeny():
    yield 20

@pytest.fixture
def nrep():
    yield 5

@pytest.fixture
def unique_parents():
    yield True

@pytest.fixture
def mateprot():
    yield TwoWayDHCross()


@pytest.fixture
def trait_mean(ntrait):
    yield numpy.zeros(ntrait)

@pytest.fixture
def trait_cov(ntrait):
    out = numpy.random.random(ntrait)
    out = numpy.outer(out, out)
    sign = numpy.random.choice([-1,1], ntrait)
    sign = numpy.outer(sign, sign)
    out = sign * out
    numpy.fill_diagonal(out, 1)
    yield out

@pytest.fixture
def embv(ndecn, trait_mean, trait_cov):
    yield numpy.random.multivariate_normal(
        mean = trait_mean,
        cov = trait_cov,
        size = ndecn
    )

@pytest.fixture
def ndecn(ntaxa):
    yield (ntaxa * (ntaxa-1)) // 2

@pytest.fixture
def decn_space_lower(ndecn):
    yield numpy.repeat(0, ndecn)

@pytest.fixture
def decn_space_upper(ndecn):
    yield numpy.repeat(1, ndecn)

@pytest.fixture
def decn_space(decn_space_lower, decn_space_upper):
    yield numpy.stack([decn_space_lower, decn_space_upper])

@pytest.fixture
def nobj():
    yield 2

@pytest.fixture
def obj_wt():
    yield numpy.array([1,-1], dtype=float)

@pytest.fixture
def obj_trans():
    yield None

@pytest.fixture
def obj_trans_kwargs():
    yield None

@pytest.fixture
def nineqcv():
    yield 0

@pytest.fixture
def ineqcv_wt():
    yield numpy.array([], dtype=float)

@pytest.fixture
def ineqcv_trans():
    yield None

@pytest.fixture
def ineqcv_trans_kwargs():
    yield None

@pytest.fixture
def neqcv():
    yield 0

@pytest.fixture
def eqcv_wt():
    yield numpy.array([], dtype=float)

@pytest.fixture
def eqcv_trans():
    yield None

@pytest.fixture
def eqcv_trans_kwargs():
    yield None

@pytest.fixture
def prob(
    embv, 
    ndecn, decn_space, decn_space_lower, decn_space_upper, 
    nobj, obj_wt, obj_trans, obj_trans_kwargs, 
    nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
    neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
):
    yield ExpectedMaximumBreedingValueBinarySelectionProblem(
        embv = embv,
        ndecn = ndecn,
        decn_space = decn_space,
        decn_space_lower = decn_space_lower,
        decn_space_upper = decn_space_upper,
        nobj = nobj,
        obj_wt = obj_wt,
        obj_trans = obj_trans,
        obj_trans_kwargs = obj_trans_kwargs,
        nineqcv = nineqcv,
        ineqcv_wt = ineqcv_wt,
        ineqcv_trans = ineqcv_trans,
        ineqcv_trans_kwargs = ineqcv_trans_kwargs,
        neqcv = neqcv,
        eqcv_wt = eqcv_wt,
        eqcv_trans = eqcv_trans,
        eqcv_trans_kwargs = eqcv_trans_kwargs
    )

########### Breeding value matrix class fixtures ###########
@pytest.fixture
def bvmat(embv):
    yield DenseBreedingValueMatrix(
        mat = embv,
        location = 0.0,
        scale = 1.0
    )

@pytest.fixture
def pgmat(ntaxa, nvrnt):
    mat = numpy.random.binomial(1, 0.5, (2, ntaxa, nvrnt)).astype("int8")
    taxa = numpy.array(["Taxa"+str(i).zfill(3) for i in range(1,ntaxa+1)], dtype=object)
    vrnt_xoprob = numpy.random.uniform(0.0, 0.1, nvrnt)
    yield DensePhasedGenotypeMatrix(
        mat = mat,
        taxa = taxa,
        vrnt_xoprob = vrnt_xoprob
    )

@pytest.fixture
def gmat(ntaxa, nvrnt):
    mat = numpy.random.binomial(2, 0.5, (ntaxa, nvrnt)).astype("int8")
    taxa = numpy.array(["Taxa"+str(i).zfill(3) for i in range(1,ntaxa+1)], dtype=object)
    yield DenseGenotypeMatrix(
        mat = mat,
        taxa = taxa
    )

@pytest.fixture
def gpmod(nvrnt, ntrait, trait_mean, trait_cov):
    beta = numpy.ones((1,ntrait), dtype=float)
    u_a = numpy.random.multivariate_normal(trait_mean, trait_cov, (nvrnt,))
    trait = numpy.array(["Trait"+str(i).zfill(2) for i in range(1,ntrait+1)], dtype=object)
    yield DenseAdditiveLinearGenomicModel(
        beta = beta,
        u_misc = None,
        u_a = u_a,
        trait = trait
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_BinaryConventionalSelectionProblem_docstring():
    assert_class_documentation(ExpectedMaximumBreedingValueBinarySelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_property_isconcrete(ExpectedMaximumBreedingValueBinarySelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == ntrait

############
### embv ###
############
def test_embv_is_concrete():
    assert_property_isconcrete(ExpectedMaximumBreedingValueBinarySelectionProblem, "embv")

def test_embv_fget(prob, ndecn, ntrait):
    assert isinstance(prob.embv, numpy.ndarray)
    assert prob.embv.shape == (ndecn,ntrait)

def test_embv_fset(prob, ndecn, ntrait):
    with not_raises(Exception):
        prob.embv = numpy.random.random((ndecn,ntrait))

def test_embv_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.embv = None
    with pytest.raises(TypeError):
        prob.embv = "string"
    with pytest.raises(TypeError):
        prob.embv = int(1)
    with pytest.raises(TypeError):
        prob.embv = float(1.0)

def test_embv_fset_ValueError(prob, ndecn, ntrait):
    with pytest.raises(ValueError):
        prob.embv = numpy.random.random(ndecn)
    with pytest.raises(ValueError):
        prob.embv = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.embv = numpy.random.random((ndecn,ndecn,ntrait,ntrait))

def test_embv_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.embv

################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test_init_is_concrete():
    assert_method_isconcrete(ExpectedMaximumBreedingValueBinarySelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_method_isconcrete(prob, "latentfn")

def test_latentfn(prob, ndecn, embv):
    x = numpy.random.binomial(1, 0.5, ndecn)
    y = (1.0 / x.sum()) * x
    a = prob.latentfn(x)
    b = -y.dot(embv)
    assert numpy.all(numpy.isclose(a,b))

################################################################################
############################## Test class methods ##############################
################################################################################
def test_from_pgmat_gpmod(
        nparent, ncross, nprogeny, nrep, unique_parents, pgmat, gpmod, mateprot,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct problem
    gebvprob = ExpectedMaximumBreedingValueBinarySelectionProblem.from_pgmat_gpmod(
        nparent, ncross, nprogeny, nrep, unique_parents, pgmat, gpmod, mateprot,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.binomial(1, 0.5, ndecn)
    y = (1.0 / x.sum()) * x
    a = gebvprob.latentfn(x)
    b = -y.dot(gebvprob.embv)
    assert numpy.all(numpy.isclose(a,b))
