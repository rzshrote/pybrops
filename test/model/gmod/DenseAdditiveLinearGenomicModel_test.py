import copy
from pathlib import Path
import numpy
import pandas
import pytest
import os
import h5py

from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import check_is_DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.test.assert_python import assert_classmethod_isconcrete, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
################# General shape parameters #################
############################################################
@pytest.fixture
def ntrait():
    yield 3

@pytest.fixture
def ntaxa():
    yield 100

@pytest.fixture
def nmarker():
    yield 100

@pytest.fixture
def nphase():
    yield 2

@pytest.fixture
def nchrom():
    yield 10

@pytest.fixture
def ngroup():
    yield 5

############################################################
###################### Genomic model #######################
############################################################
@pytest.fixture
def algmod_beta(ntrait):
    out = numpy.random.uniform(0, 10, (1,ntrait))
    yield out

@pytest.fixture
def algmod_u_misc():
    yield None

@pytest.fixture
def algmod_u_a(nmarker, ntrait):
    out = numpy.random.normal(size = (nmarker, ntrait))
    yield out

@pytest.fixture
def algmod_trait(ntrait):
    out = numpy.array(["Trait"+str(i+1) for i in range(ntrait)], dtype=object)
    yield out

@pytest.fixture
def algmod_model_name():
    yield "test_algmod"

@pytest.fixture
def algmod_params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def algmod(
        algmod_beta, 
        algmod_u_misc, 
        algmod_u_a, 
        algmod_trait, 
        algmod_model_name, 
        algmod_params
    ):
    yield DenseAdditiveLinearGenomicModel(
        beta       = algmod_beta,
        u_misc     = algmod_u_misc,
        u_a        = algmod_u_a,
        trait      = algmod_trait,
        model_name = algmod_model_name,
        hyperparams     = algmod_params
    )

############################################################
######################## Genotypes #########################
############################################################
@pytest.fixture
def pgmat_mat(nphase, ntaxa, nmarker):
    out = numpy.random.randint(0,2,(nphase,ntaxa,nmarker))
    out = out.astype("int8")
    yield out

@pytest.fixture
def pgmat_chrgrp(nchrom, nmarker):
    out = numpy.random.randint(1, nchrom+1, nmarker)
    out.sort()
    yield out

@pytest.fixture
def pgmat_phypos(nmarker):
    out = numpy.random.randint(1, 2**32-1, nmarker)
    out.sort()
    yield out

@pytest.fixture
def pgmat_taxa(ntaxa):
    out = numpy.array(["Line"+str(i+1).zfill(3) for i in range(ntaxa)], dtype=object)
    yield out

@pytest.fixture
def pgmat_taxa_grp(ngroup, ntaxa):
    out = numpy.random.randint(1, ngroup+1, ntaxa)
    out.sort()
    yield out

@pytest.fixture
def pgmat(
        pgmat_mat, 
        pgmat_chrgrp, 
        pgmat_phypos, 
        pgmat_taxa, 
        pgmat_taxa_grp
    ):
    yield DensePhasedGenotypeMatrix(
        mat = pgmat_mat,
        vrnt_chrgrp = pgmat_chrgrp,
        vrnt_phypos = pgmat_phypos,
        taxa = pgmat_taxa,
        taxa_grp = pgmat_taxa_grp
    )

############################################################
##################### Breeding values ######################
############################################################
@pytest.fixture
def mat_intercept(pgmat, algmod_beta):
    n = pgmat.ntaxa
    q = algmod_beta.shape[0]
    a = numpy.ones((n,1), dtype = "float64")
    yield a

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(DenseAdditiveLinearGenomicModel)

################################################################################
############################ Test Class Properties #############################
################################################################################

def test_u_misc_fget(algmod, algmod_u_misc):
    assert numpy.all(algmod.u_misc == algmod_u_misc)

def test_u_a_fget(algmod, algmod_u_a):
    assert numpy.all(algmod.u_a == algmod_u_a)

def test_beta_fget(algmod, algmod_beta):
    assert numpy.all(algmod.beta == algmod_beta)

def test_trait_fget(algmod, algmod_trait):
    assert numpy.all(algmod.trait == algmod_trait)

def test_algmod_model_name_fget(algmod, algmod_model_name):
    assert algmod.model_name == algmod_model_name

def test_algmod_params_fget(algmod, algmod_params):
    assert algmod.hyperparams == algmod_params

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

### __init__
def test___init___is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "__init__")

### __copy__
def test___copy___is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "__copy__")

def test___copy__(algmod):
    tmp = algmod.__copy__()
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams
    tmp = copy.copy(algmod)
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams

### __deepcopy__
def test___deepcopy___is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "__deepcopy__")

def test___deepcopy__(algmod):
    tmp = algmod.__deepcopy__()
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams
    tmp = copy.deepcopy(algmod)
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams

################################################################################
############################# Test concrete methods ############################
################################################################################

############################################################
####################### Copy methods #######################
############################################################

### copy
def test_copy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "copy")

def test_copy(algmod):
    tmp = algmod.copy()
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams

### deepcopy
def test_deepcopy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "deepcopy")

def test_deepcopy(algmod):
    tmp = algmod.deepcopy()
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams

############################################################
#################### Prediction methods ####################
############################################################

### fit_numpy
def test_fit_numpy_is_concrete():
    assert_classmethod_isconcrete(DenseAdditiveLinearGenomicModel, "fit_numpy")

def test_fit_numpy(algmod):
    with pytest.raises(AttributeError):
        algmod.fit_numpy(None, None, None)

### fit
def test_fit_is_concrete():
    assert_classmethod_isconcrete(DenseAdditiveLinearGenomicModel, "fit")

def test_fit(algmod):
    with pytest.raises(AttributeError):
        algmod.fit(None, None, None)

### predict_numpy
def test_predict_numpy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "predict_numpy")

def test_predict_numpy(algmod, mat_intercept, algmod_beta, pgmat_mat, algmod_u_a):
    geno = pgmat_mat.sum(0)
    a = algmod.predict_numpy(mat_intercept, geno)
    b = (mat_intercept @ algmod_beta) + (geno @ algmod_u_a)
    assert numpy.all(a == b)
    assert isinstance(a, numpy.ndarray)

### predict
def test_predict_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "predict")

# def test_predict(algmod, mat_intercept, algmod_beta, pgmat_mat, algmod_u_a, pgmat):
#     geno = pgmat_mat.sum(0)
#     a = algmod.predict(mat_intercept, pgmat)
#     b = (mat_intercept @ algmod_beta) + (geno @ algmod_u_a)
#     b = (b - b.mean(0)) / b.std(0)
#     assert numpy.all(a == b)
#     assert isinstance(a, BreedingValueMatrix)

### score_numpy
def test_score_numpy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "score_numpy")

def test_score_numpy(algmod, mat_intercept, pgmat_mat):
    geno = pgmat_mat.sum(0)
    y_true = algmod.predict_numpy(mat_intercept, geno)
    out = algmod.score_numpy(y_true, mat_intercept, geno)
    assert isinstance(out, numpy.ndarray)
    assert numpy.all(out == 1.0)
    assert len(out) == algmod.ntrait

### score
def test_score_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "score")

def test_score(algmod, mat_intercept, pgmat):
    y_true = algmod.predict(mat_intercept, pgmat)
    out = algmod.score(y_true, mat_intercept, pgmat)
    assert isinstance(out, numpy.ndarray)
    assert numpy.all(out == 1.0)

### gebv_numpy
def test_gebv_numpy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "gebv_numpy")

def test_gebv_numpy(algmod, pgmat_mat):
    geno = pgmat_mat.sum(0)
    out = algmod.gebv_numpy(geno)
    assert isinstance(out, numpy.ndarray)

def test_gebv_numpy_TypeError(algmod):
    with pytest.raises(TypeError):
        tmp = algmod.gebv_numpy(Z = object())

def test_gebv_numpy_ValueError(algmod):
    with pytest.raises(ValueError):
        tmp = algmod.gebv_numpy(Z = numpy.ndarray((0,)))
    with pytest.raises(ValueError):
        tmp = algmod.gebv_numpy(Z = numpy.ndarray((0,0)))

### gebv
def test_gebv_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "gebv")

def test_gebv(algmod, pgmat, pgmat_mat):
    out = algmod.gebv(pgmat_mat.sum(0))
    assert isinstance(out, BreedingValueMatrix)
    out = algmod.gebv(pgmat)
    assert isinstance(out, BreedingValueMatrix)

def test_gebv_TypeError(algmod):
    with pytest.raises(TypeError):
        tmp = algmod.gebv(Z = object())

def test_gebv_ValueError(algmod):
    with pytest.raises(ValueError):
        tmp = algmod.gebv_numpy(Z = numpy.ndarray((0,)))
    with pytest.raises(ValueError):
        tmp = algmod.gebv_numpy(Z = numpy.ndarray((0,0)))

### gegv_numpy
def test_gegv_numpy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "gegv_numpy")

def test_gegv_numpy(algmod, pgmat_mat):
    geno = pgmat_mat.sum(0)
    out = algmod.gegv_numpy(geno)
    assert isinstance(out, numpy.ndarray)

def test_gegv_numpy_TypeError(algmod):
    with pytest.raises(TypeError):
        tmp = algmod.gegv_numpy(Z = object())

def test_gegv_numpy_ValueError(algmod):
    with pytest.raises(ValueError):
        tmp = algmod.gegv_numpy(Z = numpy.ndarray((0,)))
    with pytest.raises(ValueError):
        tmp = algmod.gegv_numpy(Z = numpy.ndarray((0,0)))

### gegv
def test_gegv_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "gegv")

def test_gegv(algmod, pgmat, pgmat_mat):
    out = algmod.gegv(pgmat_mat.sum(0))
    assert isinstance(out, BreedingValueMatrix)
    out = algmod.gegv(pgmat)
    assert isinstance(out, BreedingValueMatrix)

def test_gegv_TypeError(algmod):
    with pytest.raises(TypeError):
        tmp = algmod.gegv(Z = object())

def test_gegv_ValueError(algmod):
    with pytest.raises(ValueError):
        tmp = algmod.gegv_numpy(Z = numpy.ndarray((0,)))
    with pytest.raises(ValueError):
        tmp = algmod.gegv_numpy(Z = numpy.ndarray((0,0)))

############################################################
################ Variance calculation tests ################
############################################################

### var_G_numpy
def test_var_G_numpy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "var_G_numpy")

def test_var_G_numpy(algmod, pgmat_mat):
    geno = pgmat_mat.sum(0)
    out = algmod.var_G_numpy(geno)
    assert isinstance(out, numpy.ndarray)

def test_var_G_numpy_TypeError(algmod):
    with pytest.raises(TypeError):
        tmp = algmod.var_G_numpy(Z = object())

def test_var_G_numpy_ValueError(algmod):
    with pytest.raises(ValueError):
        tmp = algmod.var_G_numpy(Z = numpy.ndarray((0,)))
    with pytest.raises(ValueError):
        tmp = algmod.var_G_numpy(Z = numpy.ndarray((0,0)))

### var_G
def test_var_G_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "var_G")

### var_A_numpy
def test_var_A_numpy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "var_A_numpy")

def test_var_A_numpy(algmod, pgmat_mat):
    geno = pgmat_mat.sum(0)
    out = algmod.var_A_numpy(geno)
    assert isinstance(out, numpy.ndarray)

def test_var_A_numpy_TypeError(algmod):
    with pytest.raises(TypeError):
        tmp = algmod.var_A_numpy(Z = object())

def test_var_A_numpy_ValueError(algmod):
    with pytest.raises(ValueError):
        tmp = algmod.var_A_numpy(Z = numpy.ndarray((0,)))
    with pytest.raises(ValueError):
        tmp = algmod.var_A_numpy(Z = numpy.ndarray((0,0)))

### var_A
def test_var_A_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "var_A")

### var_a_numpy
def test_var_a_numpy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "var_a_numpy")

### var_a
def test_var_a_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "var_a")

### bulmer_numpy
def test_bulmer_numpy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "bulmer_numpy")

### bulmer
def test_bulmer_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "bulmer")

############################################################
################## Selection limit tests ###################
############################################################

### usl_numpy
def test_usl_numpy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "usl_numpy")

### usl
def test_usl_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "usl")

def test_usl(algmod, pgmat, algmod_u_a, pgmat_mat):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = algmod_u_a.T[:,None,None,:] * pgmat_mat
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.max((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)
    a_usl = uhaplo.sum(1)
    b_usl = algmod.usl(pgmat)

    assert numpy.all(a_usl == b_usl)

### lsl_numpy
def test_lsl_numpy_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "lsl_numpy")

### lsl
def test_lsl_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "lsl")

def test_lsl(algmod, pgmat, algmod_u_a, pgmat_mat):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = algmod_u_a.T[:,None,None,:] * pgmat_mat
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.min((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)    pass

    a_lsl = uhaplo.sum(1)
    b_lsl = algmod.lsl(pgmat)

    assert numpy.all(a_lsl == b_lsl)

############################################################
################# Allele attribute tests ###################
############################################################

### facount
def test_facount_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "facount")

def test_facount_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.facount(object())
    with pytest.raises(TypeError):
        algmod.facount(pgmat, object())

def test_facount(algmod, pgmat, nmarker, ntrait):
    out = algmod.facount(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)
    dom = pgmat.mat.sum((0,1))[:,None]
    dom = numpy.broadcast_to(dom, (nmarker,ntrait))
    rec = (1 - pgmat.mat).sum((0,1))[:,None]
    rec = numpy.broadcast_to(rec, (nmarker,ntrait))
    tmp = numpy.zeros((nmarker,ntrait), int)
    mask = algmod.u_a > 0.0
    tmp[mask] = dom[mask]
    mask = algmod.u_a < 0.0
    tmp[mask] = rec[mask]
    assert numpy.all(out == tmp)

def test_facount_dtype(algmod, pgmat):
    out = algmod.facount(pgmat, float)
    assert out.dtype == float

### fafreq
def test_fafreq_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "fafreq")

def test_fafreq_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.fafreq(object())
    with pytest.raises(TypeError):
        algmod.fafreq(pgmat, object())

def test_fafreq(algmod, pgmat, nmarker, ntrait):
    out = algmod.fafreq(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)

def test_fafreq_dtype(algmod, pgmat):
    out = algmod.fafreq(pgmat, float)
    assert out.dtype == float

### faavail
def test_faavail_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "faavail")

def test_faavail_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.faavail(object())
    with pytest.raises(TypeError):
        algmod.faavail(pgmat, object())

def test_faavail(algmod, pgmat, nmarker, ntrait):
    out = algmod.faavail(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)

def test_faavail_dtype(algmod, pgmat):
    out = algmod.faavail(pgmat, float)
    assert out.dtype == float

### fafixed
def test_fafixed_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "fafixed")

def test_fafixed_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.fafixed(object())
    with pytest.raises(TypeError):
        algmod.fafixed(pgmat, object())

def test_fafixed(algmod, pgmat, nmarker, ntrait):
    out = algmod.fafixed(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)

def test_fafixed_dtype(algmod, pgmat):
    out = algmod.fafixed(pgmat, float)
    assert out.dtype == float

### fapoly
def test_fapoly_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "fapoly")

def test_fapoly_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.fapoly(object())
    with pytest.raises(TypeError):
        algmod.fapoly(pgmat, object())

def test_fapoly(algmod, pgmat, nmarker, ntrait):
    out = algmod.fapoly(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)

def test_fapoly_dtype(algmod, pgmat):
    out = algmod.fapoly(pgmat, float)
    assert out.dtype == float

### nafixed
def test_nafixed_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "nafixed")

def test_nafixed_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.nafixed(object())
    with pytest.raises(TypeError):
        algmod.nafixed(pgmat, object())

def test_nafixed(algmod, pgmat, nmarker, ntrait):
    out = algmod.nafixed(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)

def test_nafixed_dtype(algmod, pgmat):
    out = algmod.nafixed(pgmat, float)
    assert out.dtype == float

### napoly
def test_napoly_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "napoly")

def test_napoly_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.napoly(object())
    with pytest.raises(TypeError):
        algmod.napoly(pgmat, object())

def test_napoly(algmod, pgmat, nmarker, ntrait):
    out = algmod.napoly(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)

def test_napoly_dtype(algmod, pgmat):
    out = algmod.napoly(pgmat, float)
    assert out.dtype == float

### dacount
def test_dacount_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "dacount")

def test_dacount_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.dacount(object())
    with pytest.raises(TypeError):
        algmod.dacount(pgmat, object())

def test_dacount(algmod, pgmat, nmarker, ntrait):
    out = algmod.dacount(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)

def test_dacount_dtype(algmod, pgmat):
    out = algmod.dacount(pgmat, float)
    assert out.dtype == float

### dafreq
def test_dafreq_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "dafreq")

def test_dafreq_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.dafreq(object())
    with pytest.raises(TypeError):
        algmod.dafreq(pgmat, object())

def test_dafreq(algmod, pgmat, nmarker, ntrait):
    out = algmod.dafreq(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)

def test_dafreq_dtype(algmod, pgmat):
    out = algmod.dafreq(pgmat, float)
    assert out.dtype == float

### daavail
def test_daavail_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "daavail")

def test_daavail_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.daavail(object())
    with pytest.raises(TypeError):
        algmod.daavail(pgmat, object())

def test_daavail(algmod, pgmat, nmarker, ntrait):
    out = algmod.daavail(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)

def test_daavail_dtype(algmod, pgmat):
    out = algmod.daavail(pgmat, float)
    assert out.dtype == float

### dafixed
def test_dafixed_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "dafixed")

def test_dafixed_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.dafixed(object())
    with pytest.raises(TypeError):
        algmod.dafixed(pgmat, object())

def test_dafixed(algmod, pgmat, nmarker, ntrait):
    out = algmod.dafixed(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)

def test_dafixed_dtype(algmod, pgmat):
    out = algmod.dafixed(pgmat, float)
    assert out.dtype == float

### dapoly
def test_dapoly_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "dapoly")

def test_dapoly_TypeError(algmod, pgmat):
    with pytest.raises(TypeError):
        algmod.dapoly(object())
    with pytest.raises(TypeError):
        algmod.dapoly(pgmat, object())

def test_dapoly(algmod, pgmat, nmarker, ntrait):
    out = algmod.dapoly(pgmat)
    assert isinstance(out, numpy.ndarray)
    assert out.shape == (nmarker, ntrait)

def test_dapoly_dtype(algmod, pgmat):
    out = algmod.dapoly(pgmat, float)
    assert out.dtype == float

############################################################
##################### Model I/O tests ######################
############################################################

### to_pandas_dict
def test_to_pandas_dict_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "to_pandas_dict")

def test_to_pandas_dict(algmod):
    out = algmod.to_pandas_dict()
    assert isinstance(out, dict)
    assert "beta" in out
    assert "u_misc" in out
    assert "u_a" in out
    assert isinstance(out["beta"], pandas.DataFrame)
    assert isinstance(out["u_misc"], pandas.DataFrame)
    assert isinstance(out["u_a"], pandas.DataFrame)

### to_csv_dict
def test_to_csv_dict_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "to_csv_dict")

def test_to_csv_dict(algmod):
    filenames = {
        "beta": "saved_beta.csv",
        "u_misc": "saved_u_misc.csv",
        "u_a": "saved_u_a.csv",
    }
    algmod.to_csv_dict(filenames)
    assert os.path.isfile(filenames["beta"])
    assert os.path.isfile(filenames["u_misc"])
    assert os.path.isfile(filenames["u_a"])

### to_hdf5
def test_to_hdf5_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "to_hdf5")

def test_to_hdf5(algmod, shared_datadir):
    algmod.to_hdf5(shared_datadir / "algmod.hdf5")
    algmod.to_hdf5(shared_datadir / "algmod.hdf5", "prefix")

    # test whether file was created
    assert os.path.isfile(shared_datadir / "algmod.hdf5")

    algmod1 = DenseAdditiveLinearGenomicModel.from_hdf5(shared_datadir / "algmod.hdf5")
    algmod2 = DenseAdditiveLinearGenomicModel.from_hdf5(
        shared_datadir / "algmod.hdf5",
        "prefix"
    )

    # test whether data was loaded properly
    assert numpy.all(algmod.beta == algmod1.beta)
    assert numpy.all(algmod.u_misc == algmod1.u_misc)
    assert numpy.all(algmod.u_a == algmod1.u_a)
    assert numpy.all(algmod.trait == algmod1.trait)
    assert algmod.model_name == algmod1.model_name
    assert algmod.hyperparams == algmod1.hyperparams

    assert numpy.all(algmod.beta == algmod2.beta)
    assert numpy.all(algmod.u_misc == algmod2.u_misc)
    assert numpy.all(algmod.u_a == algmod2.u_a)
    assert numpy.all(algmod.trait == algmod2.trait)
    assert algmod.model_name == algmod2.model_name
    assert algmod.hyperparams == algmod2.hyperparams

### from_pandas_dict
def test_from_pandas_dict_is_concrete():
    assert_classmethod_isconcrete(DenseAdditiveLinearGenomicModel, "from_pandas_dict")

def test_from_pandas_dict(algmod):
    df_dict = algmod.to_pandas_dict()
    out = DenseAdditiveLinearGenomicModel.from_pandas_dict(df_dict)
    assert numpy.all(out.beta == algmod.beta)
    assert numpy.all(out.u_misc == algmod.u_misc)
    assert numpy.all(out.u_a == algmod.u_a)
    assert numpy.all(out.trait == algmod.trait)

### from_csv_dict
def test_from_csv_dict_is_concrete():
    assert_classmethod_isconcrete(DenseAdditiveLinearGenomicModel, "from_csv_dict")

def test_from_csv(algmod):
    filenames = {
        "beta": "saved_beta.csv",
        "u_misc": "saved_u_misc.csv",
        "u_a": "saved_u_a.csv",
    }
    algmod.to_csv_dict(filenames)
    out = DenseAdditiveLinearGenomicModel.from_csv_dict(
        filenames = filenames,
    )
    assert isinstance(out, DenseAdditiveLinearGenomicModel)

### from_hdf5
def test_from_hdf5_is_concrete():
    assert_classmethod_isconcrete(DenseAdditiveLinearGenomicModel, "from_hdf5")

def test_from_hdf5(algmod, shared_datadir):
    algmod.to_hdf5(shared_datadir / "algmod.hdf5")
    algmod.to_hdf5(shared_datadir / "algmod.hdf5", "prefix")

    # test whether file was created
    assert os.path.isfile(shared_datadir / "algmod.hdf5")

    algmod1 = DenseAdditiveLinearGenomicModel.from_hdf5(shared_datadir / "algmod.hdf5")
    algmod2 = DenseAdditiveLinearGenomicModel.from_hdf5(
        shared_datadir / "algmod.hdf5",
        "prefix"
    )

    # test whether data was loaded properly
    assert numpy.all(algmod.beta == algmod1.beta)
    assert numpy.all(algmod.u_misc == algmod1.u_misc)
    assert numpy.all(algmod.u_a == algmod1.u_a)
    assert numpy.all(algmod.trait == algmod1.trait)
    assert algmod.model_name == algmod1.model_name
    assert algmod.hyperparams == algmod1.hyperparams

    assert numpy.all(algmod.beta == algmod2.beta)
    assert numpy.all(algmod.u_misc == algmod2.u_misc)
    assert numpy.all(algmod.u_a == algmod2.u_a)
    assert numpy.all(algmod.trait == algmod2.trait)
    assert algmod.model_name == algmod2.model_name
    assert algmod.hyperparams == algmod2.hyperparams

### to_hdf5

def test_to_hdf5_is_concrete():
    assert_method_isconcrete(DenseAdditiveLinearGenomicModel, "to_hdf5")

def test_from_hdf5_str(algmod):
    fp = "tmp.h5"
    algmod.to_hdf5(fp)
    out = DenseAdditiveLinearGenomicModel.from_hdf5(fp)
    # general
    assert algmod.nexplan == out.nexplan
    assert algmod.nparam == out.nparam
    assert algmod.model_name == out.model_name
    assert algmod.hyperparams == out.hyperparams
    # beta
    assert algmod.nexplan_beta == out.nexplan_beta
    assert algmod.nparam_beta == out.nparam_beta
    assert numpy.all(algmod.beta == out.beta)
    # u
    assert algmod.nexplan_u == out.nexplan_u
    assert algmod.nparam_u == out.nparam_u
    assert numpy.all(algmod.u == out.u)
    # u_misc
    assert algmod.nexplan_u_misc == out.nexplan_u_misc
    assert algmod.nparam_u_misc == out.nparam_u_misc
    assert numpy.all(algmod.u_misc == out.u_misc)
    # u_a
    assert algmod.nexplan_u_a == out.nexplan_u_a
    assert algmod.nparam_u_a == out.nparam_u_a
    assert numpy.all(algmod.u_a == out.u_a)
    # trait
    assert algmod.ntrait == out.ntrait
    assert numpy.all(algmod.trait == out.trait)
    os.remove(fp)

def test_from_hdf5_Path(algmod):
    fp = Path("tmp.h5")
    algmod.to_hdf5(fp)
    out = DenseAdditiveLinearGenomicModel.from_hdf5(fp)
    # general
    assert algmod.nexplan == out.nexplan
    assert algmod.nparam == out.nparam
    assert algmod.model_name == out.model_name
    assert algmod.hyperparams == out.hyperparams
    # beta
    assert algmod.nexplan_beta == out.nexplan_beta
    assert algmod.nparam_beta == out.nparam_beta
    assert numpy.all(algmod.beta == out.beta)
    # u
    assert algmod.nexplan_u == out.nexplan_u
    assert algmod.nparam_u == out.nparam_u
    assert numpy.all(algmod.u == out.u)
    # u_misc
    assert algmod.nexplan_u_misc == out.nexplan_u_misc
    assert algmod.nparam_u_misc == out.nparam_u_misc
    assert numpy.all(algmod.u_misc == out.u_misc)
    # u_a
    assert algmod.nexplan_u_a == out.nexplan_u_a
    assert algmod.nparam_u_a == out.nparam_u_a
    assert numpy.all(algmod.u_a == out.u_a)
    # trait
    assert algmod.ntrait == out.ntrait
    assert numpy.all(algmod.trait == out.trait)
    os.remove(fp)

def test_from_hdf5_h5py_File(algmod):
    fp = Path("tmp.h5")
    algmod.to_hdf5(fp)
    h5file = h5py.File(fp)
    out = DenseAdditiveLinearGenomicModel.from_hdf5(h5file)
    # general
    assert algmod.nexplan == out.nexplan
    assert algmod.nparam == out.nparam
    assert algmod.model_name == out.model_name
    assert algmod.hyperparams == out.hyperparams
    # beta
    assert algmod.nexplan_beta == out.nexplan_beta
    assert algmod.nparam_beta == out.nparam_beta
    assert numpy.all(algmod.beta == out.beta)
    # u
    assert algmod.nexplan_u == out.nexplan_u
    assert algmod.nparam_u == out.nparam_u
    assert numpy.all(algmod.u == out.u)
    # u_misc
    assert algmod.nexplan_u_misc == out.nexplan_u_misc
    assert algmod.nparam_u_misc == out.nparam_u_misc
    assert numpy.all(algmod.u_misc == out.u_misc)
    # u_a
    assert algmod.nexplan_u_a == out.nexplan_u_a
    assert algmod.nparam_u_a == out.nparam_u_a
    assert numpy.all(algmod.u_a == out.u_a)
    # trait
    assert algmod.ntrait == out.ntrait
    assert numpy.all(algmod.trait == out.trait)
    h5file.close()
    os.remove(fp)

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

### from_hdf5

def test_from_hdf5_is_concrete():
    assert_classmethod_isconcrete(DenseAdditiveLinearGenomicModel, "from_hdf5")

def test_from_hdf5_str(algmod):
    fp = "tmp.h5"
    algmod.to_hdf5(fp)
    out = DenseAdditiveLinearGenomicModel.from_hdf5(fp)
    # general
    assert algmod.nexplan == out.nexplan
    assert algmod.nparam == out.nparam
    assert algmod.model_name == out.model_name
    assert algmod.hyperparams == out.hyperparams
    # beta
    assert algmod.nexplan_beta == out.nexplan_beta
    assert algmod.nparam_beta == out.nparam_beta
    assert numpy.all(algmod.beta == out.beta)
    # u
    assert algmod.nexplan_u == out.nexplan_u
    assert algmod.nparam_u == out.nparam_u
    assert numpy.all(algmod.u == out.u)
    # u_misc
    assert algmod.nexplan_u_misc == out.nexplan_u_misc
    assert algmod.nparam_u_misc == out.nparam_u_misc
    assert numpy.all(algmod.u_misc == out.u_misc)
    # u_a
    assert algmod.nexplan_u_a == out.nexplan_u_a
    assert algmod.nparam_u_a == out.nparam_u_a
    assert numpy.all(algmod.u_a == out.u_a)
    # trait
    assert algmod.ntrait == out.ntrait
    assert numpy.all(algmod.trait == out.trait)
    os.remove(fp)

def test_from_hdf5_Path(algmod):
    fp = Path("tmp.h5")
    algmod.to_hdf5(fp)
    out = DenseAdditiveLinearGenomicModel.from_hdf5(fp)
    # general
    assert algmod.nexplan == out.nexplan
    assert algmod.nparam == out.nparam
    assert algmod.model_name == out.model_name
    assert algmod.hyperparams == out.hyperparams
    # beta
    assert algmod.nexplan_beta == out.nexplan_beta
    assert algmod.nparam_beta == out.nparam_beta
    assert numpy.all(algmod.beta == out.beta)
    # u
    assert algmod.nexplan_u == out.nexplan_u
    assert algmod.nparam_u == out.nparam_u
    assert numpy.all(algmod.u == out.u)
    # u_misc
    assert algmod.nexplan_u_misc == out.nexplan_u_misc
    assert algmod.nparam_u_misc == out.nparam_u_misc
    assert numpy.all(algmod.u_misc == out.u_misc)
    # u_a
    assert algmod.nexplan_u_a == out.nexplan_u_a
    assert algmod.nparam_u_a == out.nparam_u_a
    assert numpy.all(algmod.u_a == out.u_a)
    # trait
    assert algmod.ntrait == out.ntrait
    assert numpy.all(algmod.trait == out.trait)
    os.remove(fp)

def test_from_hdf5_h5py_File(algmod):
    fp = Path("tmp.h5")
    algmod.to_hdf5(fp)
    h5file = h5py.File(fp)
    out = DenseAdditiveLinearGenomicModel.from_hdf5(h5file)
    # general
    assert algmod.nexplan == out.nexplan
    assert algmod.nparam == out.nparam
    assert algmod.model_name == out.model_name
    assert algmod.hyperparams == out.hyperparams
    # beta
    assert algmod.nexplan_beta == out.nexplan_beta
    assert algmod.nparam_beta == out.nparam_beta
    assert numpy.all(algmod.beta == out.beta)
    # u
    assert algmod.nexplan_u == out.nexplan_u
    assert algmod.nparam_u == out.nparam_u
    assert numpy.all(algmod.u == out.u)
    # u_misc
    assert algmod.nexplan_u_misc == out.nexplan_u_misc
    assert algmod.nparam_u_misc == out.nparam_u_misc
    assert numpy.all(algmod.u_misc == out.u_misc)
    # u_a
    assert algmod.nexplan_u_a == out.nexplan_u_a
    assert algmod.nparam_u_a == out.nparam_u_a
    assert numpy.all(algmod.u_a == out.u_a)
    # trait
    assert algmod.ntrait == out.ntrait
    assert numpy.all(algmod.trait == out.trait)
    h5file.close()
    os.remove(fp)

################################################################################
######################### Test class utility functions #########################
################################################################################

### check_is_DenseAdditiveLinearGenomicModel
def test_check_is_DenseAdditiveLinearGenomicModel_is_concrete():
    assert_function_isconcrete(check_is_DenseAdditiveLinearGenomicModel)

def test_check_is_DenseAdditiveLinearGenomicModel(algmod):
    with not_raises(TypeError):
        check_is_DenseAdditiveLinearGenomicModel(algmod, "algmod")
    with pytest.raises(TypeError):
        check_is_DenseAdditiveLinearGenomicModel(None, "algmod")
