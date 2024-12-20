import numpy
import pytest
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.cmat.fcty.DenseMolecularCoancestryMatrixFactory import DenseMolecularCoancestryMatrixFactory
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_property_isconcrete

from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.breed.prot.sel.prob.L2NormGenomicSelectionProblem import L2NormGenomicSubsetSelectionProblem

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
def C(ntrait,ntaxa):
    out = numpy.empty((ntrait,ntaxa,ntaxa), dtype=float)
    for i in range(ntrait):
        X = numpy.random.random((ntaxa,ntaxa))
        K = 0.5 * X.dot(X.T)
        C = numpy.linalg.cholesky(K).T
        out[i,:,:] = C
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
def decn_space_lower(ndecn):
    yield numpy.repeat(0, ndecn)

@pytest.fixture
def decn_space_upper(ndecn):
    yield numpy.repeat(ndecn-1, ndecn)

@pytest.fixture
def decn_space(ndecn):
    yield numpy.arange(ndecn)

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
    yield L2NormGenomicSubsetSelectionProblem(
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
def mkrwt(gpmod):
    yield gpmod.u_a

@pytest.fixture
def afreq(gmat, gpmod):
    yield gpmod.fafreq(gmat)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_SubsetConventionalSelectionProblem_docstring():
    assert_class_documentation(L2NormGenomicSubsetSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_property_isconcrete(L2NormGenomicSubsetSelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == ntrait

#########
### C ###
#########
def test_C_is_concrete():
    assert_property_isconcrete(L2NormGenomicSubsetSelectionProblem, "C")

def test_C_fget(prob, ntaxa, ntrait):
    assert isinstance(prob.C, numpy.ndarray)
    assert prob.C.shape == (ntrait,ntaxa,ntaxa)
    assert numpy.all(numpy.triu(prob.C) == prob.C)

def test_C_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.C = numpy.triu(numpy.random.random((ntrait,ntaxa,ntaxa)))

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
def test___init___is_concrete():
    assert_method_isconcrete(L2NormGenomicSubsetSelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_method_isconcrete(L2NormGenomicSubsetSelectionProblem, "latentfn")

def test_latentfn(prob, ndecn, C):
    x = numpy.random.binomial(1, 0.5, ndecn)
    y = (1.0 / x.sum()) * x
    x = numpy.flatnonzero(x)
    a = prob.latentfn(x)
    numpy.random.shuffle(x)
    b = prob.latentfn(x)
    c = numpy.linalg.norm(C.dot(y), ord = 2, axis = 1)
    assert numpy.all(numpy.isclose(a,b))
    assert numpy.all(numpy.isclose(a,c))

################################################################################
############################## Test class methods ##############################
################################################################################
def test_from_gmat(
        gmat, cmatfcty, mkrwt, afreq,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct problem
    l1gsprob = L2NormGenomicSubsetSelectionProblem.from_gmat(
        gmat, cmatfcty, mkrwt, afreq,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.binomial(1, 0.5, ndecn)
    y = (1.0 / x.sum()) * x
    x = numpy.flatnonzero(x)
    a = l1gsprob.latentfn(x)
    numpy.random.shuffle(x)
    b = l1gsprob.latentfn(x)
    c = numpy.linalg.norm(l1gsprob.C.dot(y), ord = 2, axis = 1)
    assert numpy.all(numpy.isclose(a,b))
    assert numpy.all(numpy.isclose(a,c))
