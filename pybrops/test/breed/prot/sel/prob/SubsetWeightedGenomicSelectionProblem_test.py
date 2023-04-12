from numbers import Real
import numpy
import pytest

from pybrops.test.assert_python import assert_concrete_property_fget, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_property

from pybrops.breed.prot.sel.prob.SubsetWeightedGenomicSelectionProblem import SubsetWeightedGenomicSelectionProblem

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
def wgebv(Z_a_mat, u_a_mat, fafreq_mat):
    yield Z_a_mat.dot(u_a_mat * numpy.power(fafreq_mat, -0.5))

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
    yield SubsetWeightedGenomicSelectionProblem(
        Z_a = Z_a_mat,
        u_a = u_a_mat,
        fafreq = fafreq_mat,
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
def test_SubsetWeightedGenomicSelectionProblem_docstring():
    assert_docstring(SubsetWeightedGenomicSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test_init_is_concrete():
    assert_concrete_method(SubsetWeightedGenomicSelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_concrete_method(prob, "latentfn")

################################################################################
########################### Test abstract properties ###########################
################################################################################
