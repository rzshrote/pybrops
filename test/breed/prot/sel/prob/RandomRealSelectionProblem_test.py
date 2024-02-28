import numpy
import pytest
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix

from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_property_isconcrete

from pybrops.breed.prot.sel.prob.RandomSelectionProblem import RandomRealSelectionProblem

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
def rbv(ntaxa, trait_mean, trait_cov):
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
    yield numpy.repeat(0.0, ntaxa)

@pytest.fixture
def decn_space_upper(ntaxa):
    yield numpy.repeat(1.0, ntaxa)

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
    rbv, 
    ndecn, decn_space, decn_space_lower, decn_space_upper, 
    nobj, obj_wt, obj_trans, obj_trans_kwargs, 
    nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
    neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
):
    yield RandomRealSelectionProblem(
        rbv = rbv,
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
def bvmat(rbv):
    yield DenseBreedingValueMatrix(
        mat = rbv,
        location = 0.0,
        scale = 1.0
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_RealConventionalSelectionProblem_docstring():
    assert_class_documentation(RandomRealSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_property_isconcrete(RandomRealSelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == ntrait

############
### rbv ###
############
def test_rbv_is_concrete():
    assert_property_isconcrete(RandomRealSelectionProblem, "rbv")

def test_rbv_fget(prob, ntaxa, ntrait):
    assert isinstance(prob.rbv, numpy.ndarray)
    assert prob.rbv.shape == (ntaxa,ntrait)

def test_rbv_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.rbv = numpy.random.random((ntaxa,ntrait))

def test_rbv_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.rbv = None
    with pytest.raises(TypeError):
        prob.rbv = "string"
    with pytest.raises(TypeError):
        prob.rbv = int(1)
    with pytest.raises(TypeError):
        prob.rbv = float(1.0)

def test_rbv_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.rbv = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.rbv = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.rbv = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_rbv_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.rbv

################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test_init_is_concrete():
    assert_method_isconcrete(RandomRealSelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_method_isconcrete(RandomRealSelectionProblem, "latentfn")

def test_latentfn(prob, ndecn, rbv):
    x = numpy.random.random(ndecn)
    y = (1.0 / x.sum()) * x
    a = prob.latentfn(x)
    b = -y.dot(rbv)
    assert numpy.all(numpy.isclose(a,b))

################################################################################
############################## Test class methods ##############################
################################################################################
def test_from_object(
        ntaxa, ntrait,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct problem
    rbvprob = RandomRealSelectionProblem.from_object(
        ntaxa, ntrait,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.random(ndecn)
    y = (1.0 / x.sum()) * x
    a = rbvprob.latentfn(x)
    b = -y.dot(rbvprob.rbv)
    assert numpy.all(numpy.isclose(a,b))
