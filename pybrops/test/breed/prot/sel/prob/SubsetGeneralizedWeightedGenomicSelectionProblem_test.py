from numbers import Real
import numpy
import pytest

from pybrops.test.assert_python import assert_concrete_property_fget, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_property

from pybrops.breed.prot.sel.prob.GeneralizedWeightedGenomicEstimatedBreedingValueSelectionProblem import GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def ntaxa():
    yield 100

@pytest.fixture
def nvrnt():
    yield 5000

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
def Z_a_mat(ntaxa, nvrnt):
    yield numpy.random.binomial(2, 0.1, size = (ntaxa,nvrnt))

@pytest.fixture
def u_a_mat(nvrnt, ntrait):
    yield numpy.random.normal(0, 1, size = (nvrnt,ntrait))

@pytest.fixture
def fafreq_mat(Z_a_mat, u_a_mat):
    # (p,)
    afreq = (0.5/Z_a_mat.shape[0]) * Z_a_mat.sum(0)
    # (p,t)
    fafreq = numpy.where(u_a_mat >= 0, afreq[:,None], 1.0-afreq[:,None])
    # prevent division by zero
    fafreq[fafreq <= 0] = 1
    yield fafreq

@pytest.fixture
def alpha():
    yield 0.612

@pytest.fixture
def gwgebv(Z_a_mat, u_a_mat, fafreq_mat, alpha):
    yield Z_a_mat.dot(u_a_mat * numpy.power(fafreq_mat, -alpha))

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
        Z_a_mat,
        u_a_mat,
        fafreq_mat,
        alpha,
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
    yield GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem(
        Z_a = Z_a_mat,
        u_a = u_a_mat,
        fafreq = fafreq_mat,
        alpha = alpha,
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
def test_SubsetGeneralizedWeightedGenomicSelectionProblem_docstring():
    assert_docstring(GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_concrete_property_fget(GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == ntrait

###########
### Z_a ###
###########
def test_Z_a_is_concrete():
    assert_concrete_property(GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem, "Z_a")

def test_Z_a_fget(prob, ntaxa, nvrnt):
    assert isinstance(prob.Z_a, numpy.ndarray)
    assert prob.Z_a.shape == (ntaxa,nvrnt)

def test_Z_a_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.Z_a = numpy.random.random((ntaxa,ntrait))

def test_Z_a_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.Z_a = None
    with pytest.raises(TypeError):
        prob.Z_a = "string"
    with pytest.raises(TypeError):
        prob.Z_a = int(1)
    with pytest.raises(TypeError):
        prob.Z_a = float(1.0)

def test_Z_a_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.Z_a = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.Z_a = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.Z_a = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_Z_a_fdel(prob):
    del prob.Z_a
    with pytest.raises(AttributeError):
        prob.Z_a

###########
### u_a ###
###########
def test_u_a_is_concrete():
    assert_concrete_property(GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem, "u_a")

def test_u_a_fget(prob, nvrnt, ntrait):
    assert isinstance(prob.u_a, numpy.ndarray)
    assert prob.u_a.shape == (nvrnt,ntrait)

def test_u_a_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.u_a = numpy.random.random((ntaxa,ntrait))

def test_u_a_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.u_a = None
    with pytest.raises(TypeError):
        prob.u_a = "string"
    with pytest.raises(TypeError):
        prob.u_a = int(1)
    with pytest.raises(TypeError):
        prob.u_a = float(1.0)

def test_u_a_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.u_a = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.u_a = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.u_a = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_u_a_fdel(prob):
    del prob.u_a
    with pytest.raises(AttributeError):
        prob.u_a

##############
### fafreq ###
##############
def test_fafreq_is_concrete():
    assert_concrete_property(GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem, "fafreq")

def test_fafreq_fget(prob, nvrnt, ntrait):
    assert isinstance(prob.fafreq, numpy.ndarray)
    assert prob.fafreq.shape == (nvrnt,ntrait)

def test_fafreq_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.fafreq = numpy.random.random((ntaxa,ntrait))

def test_fafreq_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.fafreq = None
    with pytest.raises(TypeError):
        prob.fafreq = "string"
    with pytest.raises(TypeError):
        prob.fafreq = int(1)
    with pytest.raises(TypeError):
        prob.fafreq = float(1.0)

def test_fafreq_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.fafreq = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.fafreq = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.fafreq = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_fafreq_fdel(prob):
    del prob.fafreq
    with pytest.raises(AttributeError):
        prob.fafreq

#############
### alpha ###
#############
def test_alpha_is_concrete():
    assert_concrete_property(GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem, "alpha")

def test_alpha_fget(prob):
    assert isinstance(prob.alpha, Real)

def test_alpha_fset(prob):
    with not_raises(Exception):
        prob.alpha = int(1)
    with not_raises(Exception):
        prob.alpha = float(1.0)

def test_alpha_fset_TypeError(prob, ntaxa, ntrait):
    with pytest.raises(TypeError):
        prob.alpha = None
    with pytest.raises(TypeError):
        prob.alpha = "string"
    with pytest.raises(TypeError):
        prob.alpha = complex(1.0, -1.0)
    with pytest.raises(TypeError):
        prob.alpha = numpy.random.random((ntaxa,ntrait))

def test_alpha_fset_ValueError(prob):
    with pytest.raises(ValueError):
        prob.alpha = -1.0
    with pytest.raises(ValueError):
        prob.alpha = 2.0

def test_alpha_fdel(prob):
    del prob.alpha
    with pytest.raises(AttributeError):
        prob.alpha

#############
### wgebv ###
#############
def test_wgebv_is_concrete():
    assert_concrete_property(GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem, "wgebv")

def test_wgebv_fget(prob, ntaxa, ntrait):
    assert isinstance(prob.wgebv, numpy.ndarray)
    assert prob.wgebv.shape == (ntaxa,ntrait)

def test_wgebv_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.wgebv = numpy.random.random((ntaxa,ntrait))

def test_wgebv_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.wgebv = None
    with pytest.raises(TypeError):
        prob.wgebv = "string"
    with pytest.raises(TypeError):
        prob.wgebv = int(1)
    with pytest.raises(TypeError):
        prob.wgebv = float(1.0)

def test_wgebv_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.wgebv = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.wgebv = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.wgebv = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_wgebv_fdel(prob):
    del prob.wgebv
    with pytest.raises(AttributeError):
        prob.wgebv


################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test_init_is_concrete():
    assert_concrete_method(GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_concrete_method(prob, "latentfn")

def test_latentfn(prob, ntaxa, wgebv):
    x = numpy.random.choice(ntaxa, ntaxa // 2)
    a = prob.latentfn(x)
    b = wgebv[x,:].sum(0)
    assert numpy.all(a == b)

################################################################################
########################### Test abstract properties ###########################
################################################################################
