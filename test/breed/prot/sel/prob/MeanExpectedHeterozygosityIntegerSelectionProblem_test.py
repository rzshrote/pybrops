import numpy
import pytest
from pybrops.popgen.cmat.fcty.DenseMolecularCoancestryMatrixFactory import DenseMolecularCoancestryMatrixFactory
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

from pybrops.test.assert_python import assert_concrete_property_fget, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_property

from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.breed.prot.sel.prob.MeanExpectedHeterozygositySelectionProblem import MeanExpectedHeterozygosityIntegerSelectionProblem

################################ Test fixtures #################################

##################### Shared fixtures ######################
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

################## Problem class fixtures ##################
@pytest.fixture
def C(ntaxa):
    X = numpy.random.random((ntaxa,ntaxa))
    K = 0.5 * X.dot(X.T)
    out = numpy.linalg.cholesky(K).T
    yield out

@pytest.fixture
def cmatfcty():
    yield DenseMolecularCoancestryMatrixFactory()

@pytest.fixture
def unscale():
    yield True

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
    C,
    ndecn, decn_space, decn_space_lower, decn_space_upper, 
    nobj, obj_wt, obj_trans, obj_trans_kwargs, 
    nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
    neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
):
    yield MeanExpectedHeterozygosityIntegerSelectionProblem(
        C = C,
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
def gmat(ntaxa, nvrnt):
    mat = numpy.random.binomial(2, 0.5, (ntaxa, nvrnt)).astype("int8")
    taxa = numpy.array(["Taxa"+str(i).zfill(3) for i in range(1,ntaxa+1)], dtype=object)
    yield DenseGenotypeMatrix(
        mat = mat,
        taxa = taxa
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_IntegerConventionalSelectionProblem_docstring():
    assert_docstring(MeanExpectedHeterozygosityIntegerSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_concrete_property_fget(MeanExpectedHeterozygosityIntegerSelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == 1

#########
### C ###
#########
def test_C_is_concrete():
    assert_concrete_property(MeanExpectedHeterozygosityIntegerSelectionProblem, "C")

def test_C_fget(prob, ntaxa):
    assert isinstance(prob.C, numpy.ndarray)
    assert prob.C.shape == (ntaxa,ntaxa)
    assert numpy.all(numpy.triu(prob.C) == prob.C)

def test_C_fset(prob, ntaxa):
    with not_raises(Exception):
        prob.C = numpy.triu(numpy.random.random((ntaxa,ntaxa)))

def test_C_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.C = None
    with pytest.raises(TypeError):
        prob.C = "string"
    with pytest.raises(TypeError):
        prob.C = int(1)
    with pytest.raises(TypeError):
        prob.C = float(1.0)

def test_C_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.C = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.C = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.C = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))
    with pytest.raises(ValueError):
        prob.C = numpy.random.random((ntaxa,ntaxa))

def test_C_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.C


################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test_init_is_concrete():
    assert_concrete_method(MeanExpectedHeterozygosityIntegerSelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_concrete_method(prob, "latentfn")

def test_latentfn(prob, ndecn, C):
    x = numpy.random.randint(0, ndecn, ndecn)
    y = (1.0 / x.sum()) * x
    a = prob.latentfn(x)
    b = numpy.linalg.norm(C.dot(y), ord = 2, keepdims = True) - 1.0
    assert numpy.all(numpy.isclose(a,b))

################################################################################
############################## Test class methods ##############################
################################################################################
def test_from_gmat(
        gmat, cmatfcty,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct problem
    ocsprob = MeanExpectedHeterozygosityIntegerSelectionProblem.from_gmat(
        gmat, cmatfcty,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.randint(0, ndecn, ndecn)
    y = (1.0 / x.sum()) * x
    a = ocsprob.latentfn(x)
    b = numpy.linalg.norm(ocsprob.C.dot(y), ord = 2, keepdims = True) - 1.0
    assert numpy.all(numpy.isclose(a,b))