import numpy
import pytest
from pybrops.breed.prot.sel.prob.OptimalHaploidValueSelectionProblem import OptimalHaploidValueSubsetSelectionProblem
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_property_isconcrete

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
def nparent():
    yield 2

@pytest.fixture
def nhaploblk():
    yield 10

@pytest.fixture
def unique_parents():
    yield True

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
def ohvmat(ndecn, trait_mean, trait_cov):
    yield numpy.random.multivariate_normal(
        mean = trait_mean,
        cov = trait_cov,
        size = (ndecn,)
    )

@pytest.fixture
def ndecn(ntaxa, unique_parents):
    if unique_parents:
        yield (ntaxa * (ntaxa-1)) // 2
    else:
        yield (ntaxa * (ntaxa+1)) // 2

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
    ohvmat,
    ndecn, decn_space, decn_space_lower, decn_space_upper, 
    nobj, obj_wt, obj_trans, obj_trans_kwargs, 
    nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
    neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
):
    yield OptimalHaploidValueSubsetSelectionProblem(
        ohvmat = ohvmat,
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
def pgmat(ntaxa, nvrnt):
    mat = numpy.random.binomial(1, 0.5, (2, ntaxa, nvrnt)).astype("int8")
    taxa = numpy.array(["Taxa"+str(i).zfill(3) for i in range(1,ntaxa+1)], dtype=object)
    vrnt_chrgrp = numpy.repeat([1,2], nvrnt // 2)
    vrnt_phypos = numpy.random.randint(1, 1000000000, nvrnt)
    vrnt_phypos.sort()
    vrnt_genpos = numpy.random.random(nvrnt)
    vrnt_genpos.sort()
    out = DensePhasedGenotypeMatrix(
        mat = mat,
        taxa = taxa,
        vrnt_chrgrp = vrnt_chrgrp,
        vrnt_phypos = vrnt_phypos, 
        vrnt_genpos = vrnt_genpos,
    )
    out.group_vrnt()
    yield out

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
def test_SubsetConventionalSelectionProblem_docstring():
    assert_class_documentation(OptimalHaploidValueSubsetSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_property_isconcrete(OptimalHaploidValueSubsetSelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == ntrait

############
### ohvmat ###
############
def test_ohvmat_is_concrete():
    assert_property_isconcrete(OptimalHaploidValueSubsetSelectionProblem, "ohvmat")

def test_ohvmat_fget(prob, ntrait, ndecn):
    assert isinstance(prob.ohvmat, numpy.ndarray)
    assert prob.ohvmat.shape == (ndecn,ntrait)

def test_ohvmat_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.ohvmat = numpy.random.random((ntaxa,ntrait))

def test_ohvmat_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.ohvmat = None
    with pytest.raises(TypeError):
        prob.ohvmat = "string"
    with pytest.raises(TypeError):
        prob.ohvmat = int(1)
    with pytest.raises(TypeError):
        prob.ohvmat = float(1.0)

def test_ohvmat_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.ohvmat = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.ohvmat = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.ohvmat = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_ohvmat_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.ohvmat

################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test___init___is_concrete():
    assert_method_isconcrete(OptimalHaploidValueSubsetSelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_method_isconcrete(OptimalHaploidValueSubsetSelectionProblem, "latentfn")

def test_latentfn(prob, ndecn, ohvmat):
    x = numpy.random.choice(ndecn, ndecn // 2)
    a = prob.latentfn(x)
    b = -(1.0/len(x)) * ohvmat[x,:].sum(0)
    assert numpy.all(numpy.isclose(a,b))

################################################################################
############################## Test class methods ##############################
################################################################################
def test_from_pgmat_gpmod(
        nparent, nhaploblk, unique_parents, pgmat, gpmod,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct problem
    ohvprob = OptimalHaploidValueSubsetSelectionProblem.from_pgmat_gpmod(
        nparent, nhaploblk, unique_parents, pgmat, gpmod,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.choice(ndecn, ndecn // 2)
    a = ohvprob.latentfn(x)
    b = -(1.0/len(x)) * ohvprob.ohvmat[x,:].sum(0)
    assert numpy.all(numpy.isclose(a,b))
