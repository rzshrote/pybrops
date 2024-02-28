import numpy
import pytest
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_property_isconcrete

from pybrops.breed.prot.sel.prob.GenotypeBuilderSelectionProblem import GenotypeBuilderSubsetSelectionProblem

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
def nhaploblk():
    yield 10

@pytest.fixture
def nbestfndr():
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
def haplomat(ntaxa, nhaploblk, trait_mean, trait_cov):
    yield numpy.random.multivariate_normal(
        mean = trait_mean,
        cov = trait_cov,
        size = (2,ntaxa,nhaploblk)
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
    haplomat, nbestfndr, 
    ndecn, decn_space, decn_space_lower, decn_space_upper, 
    nobj, obj_wt, obj_trans, obj_trans_kwargs, 
    nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
    neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
):
    yield GenotypeBuilderSubsetSelectionProblem(
        haplomat = haplomat,
        nbestfndr = nbestfndr,
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
    assert_class_documentation(GenotypeBuilderSubsetSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_property_isconcrete(GenotypeBuilderSubsetSelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == ntrait

############
### haplomat ###
############
def test_haplomat_is_concrete():
    assert_property_isconcrete(GenotypeBuilderSubsetSelectionProblem, "haplomat")

def test_haplomat_fget(prob, nhaploblk, ntaxa, ntrait):
    assert isinstance(prob.haplomat, numpy.ndarray)
    assert prob.haplomat.shape == (2,ntaxa,nhaploblk,ntrait)

def test_haplomat_fset(prob, nhaploblk, ntaxa, ntrait):
    with not_raises(Exception):
        prob.haplomat = numpy.random.random((2,nhaploblk,ntaxa,ntrait))

def test_haplomat_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.haplomat = None
    with pytest.raises(TypeError):
        prob.haplomat = "string"
    with pytest.raises(TypeError):
        prob.haplomat = int(1)
    with pytest.raises(TypeError):
        prob.haplomat = float(1.0)

def test_haplomat_fset_ValueError(prob, nhaploblk, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.haplomat = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.haplomat = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.haplomat = numpy.random.random((nhaploblk,ntaxa,ntrait))

def test_haplomat_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.haplomat

################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test_init_is_concrete():
    assert_method_isconcrete(GenotypeBuilderSubsetSelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_method_isconcrete(GenotypeBuilderSubsetSelectionProblem, "latentfn")

def test_latentfn(prob, ntaxa, haplomat):
    x = numpy.random.binomial(1, 0.5, ntaxa)
    x = numpy.flatnonzero(x)
    a = prob.latentfn(x)
    b = -2.0 * haplomat[:,x,:,:].max((0,1)).sum(0) # OPV upper limit
    assert numpy.all(a >= b) # GB should always be worse than OPV
    # TODO: test true GB calculations

################################################################################
############################## Test class methods ##############################
################################################################################
def test_from_object(
        ntaxa,
        pgmat, gpmod, nhaploblk, nbestfndr,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct problem
    opvprob = GenotypeBuilderSubsetSelectionProblem.from_pgmat_gpmod(
        pgmat, gpmod, nhaploblk, nbestfndr,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.binomial(1, 0.5, ntaxa)
    x = numpy.flatnonzero(x)
    a = opvprob.latentfn(x)
    b = -2.0 * opvprob.haplomat[:,x,:,:].max((0,1)).sum(0) # OPV upper limit
    assert numpy.all(a >= b) # GB should always be worse than OPV
    # TODO: test true GB calculations
