import numpy
import pytest
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_property_isconcrete

from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.breed.prot.sel.prob.L1NormGenomicSelectionProblem import L1NormGenomicBinarySelectionProblem

################################ Test fixtures #################################

##################### Shared fixtures ######################
@pytest.fixture
def ntaxa():
    yield 100

@pytest.fixture
def nvrnt():
    yield 1000

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
def V(ntrait, nvrnt, ntaxa):
    yield numpy.random.uniform(-1, 1, (ntrait,nvrnt,ntaxa))

@pytest.fixture
def ndecn(ntaxa):
    yield ntaxa

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
    V, 
    ndecn, decn_space, decn_space_lower, decn_space_upper, 
    nobj, obj_wt, obj_trans, obj_trans_kwargs, 
    nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
    neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
):
    yield L1NormGenomicBinarySelectionProblem(
        V = V,
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
def tafreq(gmat):
    yield gmat.tafreq()

@pytest.fixture
def mkrwt(gpmod):
    yield gpmod.u_a

@pytest.fixture
def tfreq(gmat, gpmod):
    yield gpmod.fafreq(gmat)



################################################################################
############################## Test class docstring ############################
################################################################################
def test_BinaryConventionalSelectionProblem_docstring():
    assert_class_documentation(L1NormGenomicBinarySelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_property_isconcrete(L1NormGenomicBinarySelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == ntrait

############
### V ###
############
def test_V_is_concrete():
    assert_property_isconcrete(L1NormGenomicBinarySelectionProblem, "V")

def test_V_fget(prob, ntaxa, ntrait, nvrnt):
    assert isinstance(prob.V, numpy.ndarray)
    assert prob.V.shape == (ntrait,nvrnt,ntaxa)

def test_V_fset(prob, ntaxa, ntrait, nvrnt):
    with not_raises(Exception):
        prob.V = numpy.random.random((ntrait,nvrnt,ntaxa))

def test_V_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.V = None
    with pytest.raises(TypeError):
        prob.V = "string"
    with pytest.raises(TypeError):
        prob.V = int(1)
    with pytest.raises(TypeError):
        prob.V = float(1.0)

def test_V_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.V = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.V = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.V = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_V_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.V

################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test___init___is_concrete():
    assert_method_isconcrete(L1NormGenomicBinarySelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_method_isconcrete(L1NormGenomicBinarySelectionProblem, "latentfn")

def test_latentfn(prob, ntaxa, V):
    x = numpy.random.binomial(1, 0.5, ntaxa)
    y = (1.0 / x.sum()) * x
    a = prob.latentfn(x)
    b = prob.latentfn(3.8*x)
    c = numpy.absolute(V.dot(y)).sum(1)
    assert numpy.all(numpy.isclose(a,b))
    assert numpy.all(numpy.isclose(a,c))

################################################################################
############################## Test class methods ##############################
################################################################################
def test_from_bvmat(
        ntaxa,
        mkrwt, tafreq, tfreq,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct problem
    l1gsprob = L1NormGenomicBinarySelectionProblem.from_numpy(
        mkrwt, tafreq, tfreq,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.binomial(1, 0.5, ntaxa)
    y = (1.0 / x.sum()) * x
    a = l1gsprob.latentfn(x)
    b = l1gsprob.latentfn(5.2*x)
    c = numpy.absolute(l1gsprob.V.dot(y)).sum(1)
    assert numpy.all(numpy.isclose(a,b))
    assert numpy.all(numpy.isclose(a,c))
