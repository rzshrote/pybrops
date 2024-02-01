from numbers import Integral
import numpy
import pytest
from pybrops.breed.prot.sel.prob.MultiObjectiveGenomicSelectionProblem import MultiObjectiveGenomicSubsetSelectionProblem
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.test.assert_python import assert_concrete_property_fget, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_property

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
    geno, ploidy, mkrwt, tfreq, 
    ndecn, decn_space, decn_space_lower, decn_space_upper, 
    nobj, obj_wt, obj_trans, obj_trans_kwargs, 
    nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
    neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
):
    yield MultiObjectiveGenomicSubsetSelectionProblem(
        geno = geno,
        ploidy = ploidy,
        mkrwt = mkrwt,
        tfreq = tfreq,
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
    vrnt_chrgrp = numpy.repeat([1,2], nvrnt // 2)
    vrnt_phypos = numpy.random.randint(1, 1000000000, nvrnt)
    vrnt_phypos.sort()
    vrnt_genpos = numpy.random.random(nvrnt)
    vrnt_genpos.sort()
    out = DenseGenotypeMatrix(
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

@pytest.fixture
def geno(gmat):
    yield gmat.mat

@pytest.fixture
def ploidy(gmat):
    yield gmat.ploidy

@pytest.fixture
def mkrwt(gpmod):
    yield numpy.absolute(gpmod.u_a)

@pytest.fixture
def tfreq(gmat, gpmod):
    yield gpmod.fafreq(gmat)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_SubsetConventionalSelectionProblem_docstring():
    assert_docstring(MultiObjectiveGenomicSubsetSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

###############
### nlatent ###
###############
def test_nlatent_is_concrete():
    assert_concrete_property_fget(MultiObjectiveGenomicSubsetSelectionProblem, "nlatent")

def test_nlatent_fget(prob, ntrait):
    assert prob.nlatent == 2 * ntrait

############
### geno ###
############
def test_geno_is_concrete():
    assert_concrete_property(MultiObjectiveGenomicSubsetSelectionProblem, "geno")

def test_geno_fget(prob, ntaxa, nvrnt):
    assert isinstance(prob.geno, numpy.ndarray)
    assert prob.geno.shape == (ntaxa,nvrnt)

def test_geno_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.geno = numpy.random.random((ntaxa,ntrait))

def test_geno_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.geno = None
    with pytest.raises(TypeError):
        prob.geno = "string"
    with pytest.raises(TypeError):
        prob.geno = int(1)
    with pytest.raises(TypeError):
        prob.geno = float(1.0)

def test_geno_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.geno = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.geno = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.geno = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_geno_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.geno

############
### ploidy ###
############
def test_ploidy_is_concrete():
    assert_concrete_property(MultiObjectiveGenomicSubsetSelectionProblem, "ploidy")

def test_ploidy_fget(prob, ntaxa, nvrnt):
    assert isinstance(prob.ploidy, Integral)

def test_ploidy_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.ploidy = int(1)

def test_ploidy_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.ploidy = None
    with pytest.raises(TypeError):
        prob.ploidy = object()
    with pytest.raises(TypeError):
        prob.ploidy = "string"
    with pytest.raises(TypeError):
        prob.ploidy = float(1.0)

def test_ploidy_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.ploidy = int(0)
    with pytest.raises(ValueError):
        prob.ploidy = int(-1)

def test_ploidy_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.ploidy

############
### mkrwt ###
############
def test_mkrwt_is_concrete():
    assert_concrete_property(MultiObjectiveGenomicSubsetSelectionProblem, "mkrwt")

def test_mkrwt_fget(prob, ntrait, nvrnt):
    assert isinstance(prob.mkrwt, numpy.ndarray)
    assert prob.mkrwt.shape == (nvrnt, ntrait)

def test_mkrwt_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.mkrwt = numpy.random.random((ntaxa,ntrait))

def test_mkrwt_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.mkrwt = None
    with pytest.raises(TypeError):
        prob.mkrwt = "string"
    with pytest.raises(TypeError):
        prob.mkrwt = int(1)
    with pytest.raises(TypeError):
        prob.mkrwt = float(1.0)

def test_mkrwt_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.mkrwt = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.mkrwt = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.mkrwt = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_mkrwt_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.mkrwt

############
### tfreq ###
############
def test_tfreq_is_concrete():
    assert_concrete_property(MultiObjectiveGenomicSubsetSelectionProblem, "tfreq")

def test_tfreq_fget(prob, ntrait, nvrnt):
    assert isinstance(prob.tfreq, numpy.ndarray)
    assert prob.tfreq.shape == (nvrnt, ntrait)

def test_tfreq_fset(prob, ntaxa, ntrait):
    with not_raises(Exception):
        prob.tfreq = numpy.random.random((ntaxa,ntrait))

def test_tfreq_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.tfreq = None
    with pytest.raises(TypeError):
        prob.tfreq = "string"
    with pytest.raises(TypeError):
        prob.tfreq = int(1)
    with pytest.raises(TypeError):
        prob.tfreq = float(1.0)

def test_tfreq_fset_ValueError(prob, ntaxa, ntrait):
    with pytest.raises(ValueError):
        prob.tfreq = numpy.random.random(ntaxa)
    with pytest.raises(ValueError):
        prob.tfreq = numpy.random.random(ntrait)
    with pytest.raises(ValueError):
        prob.tfreq = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

def test_tfreq_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.tfreq

################################################################################
############################# Test concrete methods ############################
################################################################################

################
### __init__ ###
################
def test_init_is_concrete():
    assert_concrete_method(MultiObjectiveGenomicSubsetSelectionProblem, "__init__")

################
### latentfn ###
################
def test_latentfn_is_concrete(prob):
    assert_concrete_method(prob, "latentfn")

def test_latentfn(prob, ntaxa, geno):
    x = numpy.random.choice(ntaxa, ntaxa // 2)
    a = prob.latentfn(x)
    numpy.random.shuffle(x)
    b = prob.latentfn(x)
    # b = -(1.0/len(x)) * geno[x,:].sum(0)
    assert numpy.all(numpy.isclose(a,b))

################################################################################
############################## Test class methods ##############################
################################################################################
def test_from_gmat_gpmod(
        ntaxa,
        gmat, mkrwt, tfreq, gpmod,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    ):
    # construct problem
    ohvprob = MultiObjectiveGenomicSubsetSelectionProblem.from_gmat_gpmod(
        gmat, mkrwt, tfreq, gpmod,
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, obj_trans, obj_trans_kwargs, 
        nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
        neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs
    )
    # test problem calculations
    x = numpy.random.choice(ntaxa, ntaxa // 2)
    a = ohvprob.latentfn(x)
    numpy.random.shuffle(x)
    b = ohvprob.latentfn(x)
    # c = -(1.0/len(x)) * ohvprob.geno[x,:].sum(0)
    assert numpy.all(numpy.isclose(a,b))
