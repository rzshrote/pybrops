import numpy
import pytest

from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_property_isconcrete
from pybrops.breed.prot.sel.prob.GeneralizedWeightedGenomicEstimatedBreedingValueSelectionProblem import GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def ntaxa():
    yield 100

@pytest.fixture
def ntrait():
    yield 2

@pytest.fixture
def nvrnt():
    yield 1000

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
def gwgebv(ntaxa, trait_mean, trait_cov):
    yield numpy.random.multivariate_normal(
        mean = trait_mean,
        cov = trait_cov,
        size = ntaxa
    )

@pytest.fixture
def ndecn(ntaxa):
    yield ntaxa

@pytest.fixture
def decn_space_lower(ntaxa):
    yield numpy.repeat(0, ntaxa)

@pytest.fixture
def decn_space_upper(ntaxa):
    yield numpy.repeat(1, ntaxa)

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
    gwgebv, 
    ndecn, decn_space, decn_space_lower, decn_space_upper, 
    nobj, obj_wt, obj_trans, obj_trans_kwargs, 
    nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
    neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
):
    yield GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem(
        gwgebv = gwgebv,
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
def bvmat(gwgebv):
    yield DenseBreedingValueMatrix(
        mat = gwgebv,
        location = 0.0,
        scale = 1.0
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
    assert_class_documentation(GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_property_isconcrete(GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == ntrait

############
### gwgebv ###
############
def test_gwgebv_is_concrete():
    assert_property_isconcrete(GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem, "gwgebv")

def test_gwgebv_fget(prob, ntaxa, ntrait):
    assert isinstance(prob.gwgebv, numpy.ndarray)
    assert prob.gwgebv.shape == (ntaxa,ntrait)

def test_gwgebv_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.gwgebv = numpy.random.random((ntaxa,ntrait))

def test_gwgebv_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.gwgebv = None
    with pytest.raises(TypeError):
        prob.gwgebv = "string"
    with pytest.raises(TypeError):
        prob.gwgebv = int(1)
    with pytest.raises(TypeError):
        prob.gwgebv = float(1.0)

def test_gwgebv_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.gwgebv = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.gwgebv = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.gwgebv = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_gwgebv_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.gwgebv

################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test_init_is_concrete():
    assert_method_isconcrete(GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_method_isconcrete(GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem, "latentfn")

def test_latentfn(prob, ntaxa, gwgebv):
    x = numpy.random.binomial(1, 0.5, ntaxa)
    y = (1.0 / x.sum()) * x
    a = prob.latentfn(x)
    b = -y.dot(gwgebv)
    assert numpy.all(numpy.isclose(a,b))

################################################################################
############################## Test class methods ##############################
################################################################################
def test_from_numpy(
        ntaxa, nvrnt, ntrait, trait_mean, trait_cov,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct true values
    Z_a = numpy.random.binomial(2, 0.5, (ntaxa, nvrnt))
    u_a = numpy.random.multivariate_normal(trait_mean, trait_cov, (nvrnt,))
    fafreq = numpy.random.random((nvrnt,ntrait))
    alpha = 0.618
    gwgebv = Z_a.dot(u_a * numpy.power(fafreq, -alpha))
    # construct problem
    gebvprob = GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem.from_numpy(
        Z_a, u_a, fafreq, alpha,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.binomial(1, 0.5, ntaxa)
    y = (1.0 / x.sum()) * x
    a = gebvprob.latentfn(x)
    b = -y.dot(gwgebv)
    assert numpy.all(numpy.isclose(a,b))

def test_from_gmat_algpmod(
        ntaxa,
        gmat, gpmod, 
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct true values
    Z_a = gmat.mat
    u_a = gpmod.u_a
    fafreq = gpmod.fafreq(gmat)
    alpha = 0.618
    gwgebv = Z_a.dot(u_a * numpy.power(fafreq, -alpha))
    # construct problem
    gebvprob = GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem.from_gmat_algpmod(
        gmat, gpmod, alpha,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.binomial(1, 0.5, ntaxa)
    y = (1.0 / x.sum()) * x
    a = gebvprob.latentfn(x)
    b = -y.dot(gwgebv)
    assert numpy.all(numpy.isclose(a,b))
