import numpy
import pytest

from pybrops.test.assert_python import assert_concrete_property_fget, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_property

from pybrops.breed.prot.sel.prob.EstimatedBreedingValueSelectionProblem import EstimatedBreedingValueSubsetSelectionProblem

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
def ebv(ntaxa, trait_mean, trait_cov):
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
    ebv,
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
    yield EstimatedBreedingValueSubsetSelectionProblem(
        ebv = ebv,
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

################################################################################
############################## Test class docstring ############################
################################################################################
def test_SubsetConventionalSelectionProblem_docstring():
    assert_docstring(EstimatedBreedingValueSubsetSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_concrete_property_fget(EstimatedBreedingValueSubsetSelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == ntrait

############
### ebv ###
############
def test_ebv_is_concrete():
    assert_concrete_property(EstimatedBreedingValueSubsetSelectionProblem, "ebv")

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
    assert_concrete_method(EstimatedBreedingValueSubsetSelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_concrete_method(prob, "latentfn")

def test_latentfn(prob, ntaxa, ebv):
    x = numpy.random.choice(ntaxa, ntaxa // 2)
    a = prob.latentfn(x)
    b = -(1.0/len(x)) * ebv[x,:].sum(0)
    assert numpy.all(numpy.isclose(a,b))

################################################################################
########################### Test abstract properties ###########################
################################################################################
