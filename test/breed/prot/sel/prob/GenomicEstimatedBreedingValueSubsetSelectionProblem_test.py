import numpy
import pytest

from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_property_isconcrete
from pybrops.breed.prot.sel.prob.GenomicEstimatedBreedingValueSelectionProblem import GenomicEstimatedBreedingValueSubsetSelectionProblem

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
def gebv(ntaxa, trait_mean, trait_cov):
    yield numpy.random.multivariate_normal(
        mean = trait_mean,
        cov = trait_cov,
        size = ntaxa
    )

@pytest.fixture
def ndecn():
    yield 4

@pytest.fixture
def decn_space_lower():
    yield numpy.array([1,2,3,4], dtype=float)

@pytest.fixture
def decn_space_upper():
    yield numpy.array([5,6,7,8], dtype=float)

@pytest.fixture
def decn_space(decn_space_lower, decn_space_upper):
    yield numpy.concatenate([decn_space_lower, decn_space_upper])

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
    gebv,
    ndecn,
    decn_space,
    decn_space_lower,
    decn_space_upper,
    nobj,
    obj_wt,
    obj_trans,
    obj_trans_kwargs,
    nineqcv,
    ineqcv_wt,
    ineqcv_trans,
    ineqcv_trans_kwargs,
    neqcv,
    eqcv_wt,
    eqcv_trans,
    eqcv_trans_kwargs
):
    yield GenomicEstimatedBreedingValueSubsetSelectionProblem(
        gebv = gebv,
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
def bvmat(gebv):
    yield DenseBreedingValueMatrix(
        mat = gebv,
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

@pytest.fixture
def unscale():
    yield False

################################################################################
############################## Test class docstring ############################
################################################################################
def test_SubsetConventionalSelectionProblem_docstring():
    assert_class_documentation(GenomicEstimatedBreedingValueSubsetSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_property_isconcrete(GenomicEstimatedBreedingValueSubsetSelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == ntrait

############
### gebv ###
############
def test_gebv_is_concrete():
    assert_property_isconcrete(GenomicEstimatedBreedingValueSubsetSelectionProblem, "gebv")

def test_gebv_fget(prob, ntaxa, ntrait):
    assert isinstance(prob.gebv, numpy.ndarray)
    assert prob.gebv.shape == (ntaxa,ntrait)

def test_gebv_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.gebv = numpy.random.random((ntaxa,ntrait))

def test_gebv_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.gebv = None
    with pytest.raises(TypeError):
        prob.gebv = "string"
    with pytest.raises(TypeError):
        prob.gebv = int(1)
    with pytest.raises(TypeError):
        prob.gebv = float(1.0)

def test_gebv_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.gebv = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.gebv = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.gebv = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_gebv_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.gebv

################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test___init___is_concrete():
    assert_method_isconcrete(GenomicEstimatedBreedingValueSubsetSelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_method_isconcrete(GenomicEstimatedBreedingValueSubsetSelectionProblem, "latentfn")

def test_latentfn(prob, ntaxa, gebv):
    x = numpy.random.choice(ntaxa, ntaxa // 2)
    a = prob.latentfn(x)
    b = -(1.0/len(x)) * gebv[x,:].sum(0)
    assert numpy.all(numpy.isclose(a,b))

################################################################################
############################## Test class methods ##############################
################################################################################
def test_from_bvmat(
        ntaxa, gebv,
        bvmat, unscale,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct problem
    gebvprob = GenomicEstimatedBreedingValueSubsetSelectionProblem.from_bvmat(
        bvmat, unscale,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.choice(ntaxa, ntaxa // 2)
    a = gebvprob.latentfn(x)
    b = -(1.0/len(x)) * gebv[x,:].sum(0)
    assert numpy.all(numpy.isclose(a,b))

def test_from_gmat_gpmod(
        ntaxa,
        gmat, gpmod, unscale,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct problem
    gebvprob = GenomicEstimatedBreedingValueSubsetSelectionProblem.from_gmat_gpmod(
        gmat, gpmod, unscale,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # calculate GEBVs
    gebv = gpmod.gebv(gmat).mat
    # test problem calculations
    x = numpy.random.choice(ntaxa, ntaxa // 2)
    a = gebvprob.latentfn(x)
    b = -(1.0/len(x)) * gebv[x,:].sum(0)
    assert numpy.all(numpy.isclose(a,b))
