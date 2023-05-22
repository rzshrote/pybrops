import numpy
import pytest

from pybrops.test.assert_python import assert_concrete_property_fget, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_property

from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.breed.prot.sel.prob.FamilyEstimatedBreedingValueSelectionProblem import FamilyEstimatedBreedingValueBinarySelectionProblem

################################ Test fixtures #################################

##################### Shared fixtures ######################
@pytest.fixture
def ntaxa():
    yield 100

@pytest.fixture
def ntrait():
    yield 2

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

################## Problem class fixtures ##################
@pytest.fixture
def ebv(ntaxa, trait_mean, trait_cov):
    yield numpy.random.multivariate_normal(
        mean = trait_mean,
        cov = trait_cov,
        size = ntaxa
    )

@pytest.fixture
def familyid(ntaxa):
    yield numpy.repeat(numpy.arange(ntaxa // 20), 20)

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
    ebv, familyid,
    ndecn, decn_space, decn_space_lower, decn_space_upper, 
    nobj, obj_wt, obj_trans, obj_trans_kwargs, 
    nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
    neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
):
    yield FamilyEstimatedBreedingValueBinarySelectionProblem(
        ebv = ebv,
        familyid = familyid,
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
def bvmat(ebv, familyid):
    yield DenseBreedingValueMatrix(
        mat = ebv,
        location = 0.0,
        scale = 1.0,
        taxa_grp = familyid
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_BinaryConventionalSelectionProblem_docstring():
    assert_docstring(FamilyEstimatedBreedingValueBinarySelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_concrete_property_fget(FamilyEstimatedBreedingValueBinarySelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait, familyid):
    assert prob.nlatent == (ntrait + len(numpy.unique(familyid)))

############
### ebv ###
############
def test_ebv_is_concrete():
    assert_concrete_property(FamilyEstimatedBreedingValueBinarySelectionProblem, "ebv")

def test_ebv_fget(prob, ntaxa, ntrait):
    assert isinstance(prob.ebv, numpy.ndarray)
    assert prob.ebv.shape == (ntaxa,ntrait)

def test_ebv_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.ebv = numpy.random.random((ntaxa,ntrait))

def test_ebv_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.ebv = None
    with pytest.raises(TypeError):
        prob.ebv = "string"
    with pytest.raises(TypeError):
        prob.ebv = int(1)
    with pytest.raises(TypeError):
        prob.ebv = float(1.0)

def test_ebv_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.ebv = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.ebv = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.ebv = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_ebv_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.ebv

################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test_init_is_concrete():
    assert_concrete_method(FamilyEstimatedBreedingValueBinarySelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_concrete_method(prob, "latentfn")

def test_latentfn(prob, ntaxa, ebv, familyid):
    x = numpy.random.binomial(1, 0.5, ntaxa)
    x = (1.0 / x.sum()) * x
    a = prob.latentfn(x)
    b = -x.dot(ebv)
    c = -numpy.bincount(prob.familyix, x)
    d = numpy.concatenate([b,c])
    assert numpy.all(numpy.isclose(a,d))

################################################################################
############################## Test class methods ##############################
################################################################################
def test_from_bvmat(
        ntaxa, ebv,
        bvmat, 
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct problem
    ebvprob = FamilyEstimatedBreedingValueBinarySelectionProblem.from_bvmat(
        bvmat, 
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.binomial(1, 0.5, ntaxa)
    x = (1.0 / x.sum()) * x
    a = ebvprob.latentfn(x)
    b = -x.dot(ebv)
    c = -numpy.bincount(ebvprob.familyix, x)
    d = numpy.concatenate([b,c])
    assert numpy.all(numpy.isclose(a,d))
