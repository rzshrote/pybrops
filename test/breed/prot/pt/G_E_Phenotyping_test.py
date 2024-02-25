from numbers import Integral
import numpy
from numpy.random import Generator, RandomState, PCG64
import pandas
import pytest
import copy

from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete

from pybrops.breed.prot.pt.G_E_Phenotyping import G_E_Phenotyping
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
######################## Genotypes #########################
############################################################
@pytest.fixture
def mat_int8():
    yield numpy.int8([
       [[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]],

       [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0]]
    ])

@pytest.fixture
def mat_chrgrp():
    yield numpy.int64([
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2
    ])

@pytest.fixture
def mat_phypos():
    yield numpy.int64([
         1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
        11, 12, 13, 14, 15, 16, 17, 18, 19, 20
    ])

@pytest.fixture
def mat_genpos():
    yield numpy.float64([
        0.13, 0.32, 0.53, 0.54, 0.55, 0.61, 0.63, 0.7 , 0.75, 0.96,
        0.14, 0.16, 0.26, 0.31, 0.31, 0.68, 0.7 , 0.74, 0.75, 0.91
    ])

@pytest.fixture
def mat_taxa():
    yield numpy.array([
        'Line01', 'Line02', 'Line03', 'Line04', 'Line05',
        'Line06', 'Line07', 'Line08', 'Line09', 'Line10',
        'Line11', 'Line12', 'Line13', 'Line14', 'Line15',
        'Line16', 'Line17', 'Line18', 'Line19', 'Line20'
    ], dtype = object)

@pytest.fixture
def mat_taxa_grp():
    yield numpy.int64([
        1, 1, 1, 1, 1,
        2, 2, 2, 2, 2,
        3, 3, 3, 3, 3,
        4, 4, 4, 4, 4
    ])

@pytest.fixture
def dpgmat(mat_int8, mat_chrgrp, mat_phypos, mat_genpos, mat_taxa, mat_taxa_grp):
    out = DensePhasedGenotypeMatrix(
        mat = mat_int8,
        vrnt_chrgrp = mat_chrgrp,
        vrnt_phypos = mat_phypos,
        vrnt_genpos = mat_genpos,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp
    )
    out.group_vrnt()
    yield out

############################################################
###################### Genomic model #######################
############################################################
@pytest.fixture
def ntrait():
    yield 2

@pytest.fixture
def mat_beta():
    yield numpy.float64([
        [25.6, 13.4]
    ])
    # yield numpy.float64([[1.4, 2.5, 7.2]])

@pytest.fixture
def mat_u_misc():
    yield None

@pytest.fixture
def mat_u_a():
    yield numpy.float64([
        [ 1.25, -0.68],
        [-0.02, -1.09],
        [ 0.21, -0.5 ],
        [-2.84,  0.64],
        [-1.37, -0.81],
        [-2.06,  2.22],
        [ 1.52, -0.21],
        [-0.23, -1.78],
        [ 1.04, -0.55],
        [-0.77, -1.4 ],
        [-0.44,  0.89],
        [ 0.12, -0.87],
        [-0.44, -0.55],
        [ 1.36,  0.73],
        [ 1.04,  1.22],
        [-0.05,  0.82],
        [ 0.93,  0.73],
        [-0.89,  1.21],
        [ 0.05, -1.19],
        [-1.27, -2.  ]
    ])

@pytest.fixture
def trait():
    yield numpy.array(["protein", "yield"], dtype = object)

@pytest.fixture
def model_name():
    yield "test_dalgmod"

@pytest.fixture
def hyperparams():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def gpmod(mat_beta, mat_u_misc, mat_u_a, trait, model_name, hyperparams):
    yield DenseAdditiveLinearGenomicModel(
        beta = mat_beta,
        u_misc = mat_u_misc,
        u_a = mat_u_a,
        trait = trait,
        model_name = model_name,
        hyperparams = hyperparams
    )

############################################################
################## Breeding values model ###################
############################################################
@pytest.fixture
def bvmat(gpmod, dpgmat):
    yield gpmod.gebv(dpgmat)

############################################################
##################### G_E_Phenotyping ######################
############################################################
@pytest.fixture
def nenv():
    yield 2

@pytest.fixture
def nrep():
    yield 2

@pytest.fixture
def var_env():
    yield 1

@pytest.fixture
def var_rep():
    yield 0.25

@pytest.fixture
def var_err():
    yield 0.5

@pytest.fixture
def ptprot(gpmod, nenv, nrep, var_env, var_rep, var_err):
    yield G_E_Phenotyping(
        gpmod = gpmod,
        nenv = nenv,
        nrep = nrep,
        var_env = var_env,
        var_rep = var_rep,
        var_err = var_err
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(G_E_Phenotyping)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(G_E_Phenotyping, "__init__")

def test_phenotype_is_concrete():
    assert_method_isconcrete(G_E_Phenotyping, "phenotype")

def test_set_h2_is_concrete():
    assert_method_isconcrete(G_E_Phenotyping, "set_h2")

def test_set_H2_is_concrete():
    assert_method_isconcrete(G_E_Phenotyping, "set_H2")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

### gpmod ###
def test_gpmod_is_conrete():
    assert_property_isconcrete(G_E_Phenotyping, "gpmod")

def test_gpmod_fget(ptprot, gpmod):
    # TODO: assert equality of models
    assert id(ptprot.gpmod) == id(gpmod)

def test_gpmod_fset(ptprot, gpmod):
    a = copy.copy(gpmod)
    ptprot.gpmod = a
    assert id(ptprot.gpmod) == id(a)

def test_gpmod_fdel(ptprot):
    with pytest.raises(AttributeError):
        del ptprot.gpmod

### nenv ###
def test_nenv_is_concrete():
    assert_property_isconcrete(G_E_Phenotyping, "nenv")

def test_nenv_fget(ptprot, nenv):
    assert ptprot.nenv == nenv
    assert isinstance(ptprot.nenv, Integral)

def test_nenv_fset(ptprot, nenv):
    ptprot.nenv = nenv
    assert ptprot.nenv == nenv

def test_nenv_fset_TypeError(ptprot, nenv):
    with pytest.raises(TypeError):
        ptprot.nenv = object()
    with pytest.raises(TypeError):
        ptprot.nenv = float(1)
    with pytest.raises(TypeError):
        ptprot.nenv = numpy.array([nenv])
    with not_raises(TypeError):
        ptprot.nenv = int(nenv)
    with not_raises(TypeError):
        ptprot.nenv = numpy.int8(nenv)
    with not_raises(TypeError):
        ptprot.nenv = numpy.int16(nenv)
    with not_raises(TypeError):
        ptprot.nenv = numpy.int32(nenv)
    with not_raises(TypeError):
        ptprot.nenv = numpy.int64(nenv)

def test_nenv_fset_ValueError(ptprot, nenv):
    with pytest.raises(ValueError):
        ptprot.nenv = int(0)
    with pytest.raises(ValueError):
        ptprot.nenv = int(-1)
    with not_raises(ValueError):
        ptprot.nenv = nenv

def test_nenv_fdel(ptprot):
    with pytest.raises(AttributeError):
        del ptprot.nenv

### nrep ###
def test_nrep_is_concrete():
    assert_property_isconcrete(G_E_Phenotyping, "nrep")

def test_nrep_fget(ptprot, nrep):
    assert numpy.all(ptprot.nrep == nrep)
    assert isinstance(ptprot.nrep, numpy.ndarray)

def test_nrep_fset(ptprot, nrep):
    ptprot.nrep = nrep
    assert numpy.all(ptprot.nrep == nrep)

def test_nrep_fset_TypeError(ptprot, nenv, nrep):
    with pytest.raises(TypeError):
        ptprot.nrep = object()
    with pytest.raises(TypeError):
        ptprot.nrep = float(1)
    with not_raises(TypeError):
        ptprot.nrep = int(nrep)
    with not_raises(TypeError):
        ptprot.nrep = numpy.int8(nrep)
    with not_raises(TypeError):
        ptprot.nrep = numpy.int16(nrep)
    with not_raises(TypeError):
        ptprot.nrep = numpy.int32(nrep)
    with not_raises(TypeError):
        ptprot.nrep = numpy.int64(nrep)
    with not_raises(TypeError):
        ptprot.nrep = numpy.repeat(nrep, nenv)

def test_nrep_fset_ValueError(ptprot, nenv, nrep):
    with pytest.raises(ValueError):
        ptprot.nrep = int(0)
    with pytest.raises(ValueError):
        ptprot.nrep = int(-1)
    with pytest.raises(ValueError):
        ptprot.nrep = numpy.repeat(nrep, nenv-1)
    with pytest.raises(ValueError):
        ptprot.nrep = numpy.repeat(nrep, nenv+1)
    with not_raises(ValueError):
        ptprot.nrep = nrep

def test_nrep_fdel(ptprot):
    with pytest.raises(AttributeError):
        del ptprot.nrep

### var_env ###
def test_var_env_is_concrete():
    assert_property_isconcrete(G_E_Phenotyping, "var_env")

def test_var_env_fget(ptprot, var_env):
    assert numpy.all(ptprot.var_env == var_env)
    assert isinstance(ptprot.var_env, numpy.ndarray)

def test_var_env_fset(ptprot, var_env):
    ptprot.var_env = var_env
    assert numpy.all(ptprot.var_env == var_env)

def test_var_env_fset_TypeError(ptprot, ntrait, var_env):
    with pytest.raises(TypeError):
        ptprot.var_env = object()
    with not_raises(TypeError):
        ptprot.var_env = int(1)
    with not_raises(TypeError):
        ptprot.var_env = float(var_env)
    with not_raises(TypeError):
        ptprot.var_env = numpy.int8(var_env)
    with not_raises(TypeError):
        ptprot.var_env = numpy.int16(var_env)
    with not_raises(TypeError):
        ptprot.var_env = numpy.int32(var_env)
    with not_raises(TypeError):
        ptprot.var_env = numpy.int64(var_env)
    with not_raises(TypeError):
        ptprot.var_env = numpy.repeat(var_env, ntrait)

def test_var_env_fset_ValueError(ptprot, ntrait, var_env):
    with pytest.raises(ValueError):
        ptprot.var_env = float(-1)
    with pytest.raises(ValueError):
        ptprot.var_env = numpy.repeat(var_env, ntrait-1)
    with pytest.raises(ValueError):
        ptprot.var_env = numpy.repeat(var_env, ntrait+1)
    with not_raises(ValueError):
        ptprot.var_env = float(0)
    with not_raises(ValueError):
        ptprot.var_env = var_env

def test_var_env_fdel(ptprot):
    with pytest.raises(AttributeError):
        del ptprot.var_env

### var_rep ###
def test_var_rep_is_concrete():
    assert_property_isconcrete(G_E_Phenotyping, "var_rep")

def test_var_rep_fget(ptprot, var_rep):
    assert numpy.all(ptprot.var_rep == var_rep)
    assert isinstance(ptprot.var_rep, numpy.ndarray)

def test_var_rep_fset(ptprot, var_rep):
    ptprot.var_rep = var_rep
    assert numpy.all(ptprot.var_rep == var_rep)

def test_var_rep_fset_TypeError(ptprot, ntrait, var_rep):
    with pytest.raises(TypeError):
        ptprot.var_rep = object()
    with not_raises(TypeError):
        ptprot.var_rep = int(1)
    with not_raises(TypeError):
        ptprot.var_rep = float(var_rep)
    with not_raises(TypeError):
        ptprot.var_rep = numpy.int8(var_rep)
    with not_raises(TypeError):
        ptprot.var_rep = numpy.int16(var_rep)
    with not_raises(TypeError):
        ptprot.var_rep = numpy.int32(var_rep)
    with not_raises(TypeError):
        ptprot.var_rep = numpy.int64(var_rep)
    with not_raises(TypeError):
        ptprot.var_rep = numpy.repeat(var_rep, ntrait)

def test_var_rep_fset_ValueError(ptprot, ntrait, var_rep):
    with pytest.raises(ValueError):
        ptprot.var_rep = float(-1)
    with pytest.raises(ValueError):
        ptprot.var_rep = numpy.repeat(var_rep, ntrait-1)
    with pytest.raises(ValueError):
        ptprot.var_rep = numpy.repeat(var_rep, ntrait+1)
    with not_raises(ValueError):
        ptprot.var_rep = float(0)
    with not_raises(ValueError):
        ptprot.var_rep = var_rep

def test_var_rep_fdel(ptprot):
    with pytest.raises(AttributeError):
        del ptprot.var_rep

### var_err ###
def test_var_err_is_concrete():
    assert_property_isconcrete(G_E_Phenotyping, "var_err")

def test_var_err_fget(ptprot, var_err):
    assert numpy.all(ptprot.var_err == var_err)
    assert isinstance(ptprot.var_err, numpy.ndarray)

def test_var_err_fset(ptprot, var_err):
    ptprot.var_err = var_err
    assert numpy.all(ptprot.var_err == var_err)

def test_var_err_fset_TypeError(ptprot, ntrait, var_err):
    with pytest.raises(TypeError):
        ptprot.var_err = object()
    with not_raises(TypeError):
        ptprot.var_err = int(1)
    with not_raises(TypeError):
        ptprot.var_err = float(var_err)
    with not_raises(TypeError):
        ptprot.var_err = numpy.int8(var_err)
    with not_raises(TypeError):
        ptprot.var_err = numpy.int16(var_err)
    with not_raises(TypeError):
        ptprot.var_err = numpy.int32(var_err)
    with not_raises(TypeError):
        ptprot.var_err = numpy.int64(var_err)
    with not_raises(TypeError):
        ptprot.var_err = numpy.repeat(var_err, ntrait)

def test_var_err_fset_ValueError(ptprot, ntrait, var_err):
    with pytest.raises(ValueError):
        ptprot.var_err = float(-1)
    with pytest.raises(ValueError):
        ptprot.var_err = numpy.repeat(var_err, ntrait-1)
    with pytest.raises(ValueError):
        ptprot.var_err = numpy.repeat(var_err, ntrait+1)
    with not_raises(ValueError):
        ptprot.var_err = float(0)
    with not_raises(ValueError):
        ptprot.var_err = var_err

def test_var_err_fdel(ptprot):
    with pytest.raises(AttributeError):
        del ptprot.var_err

### rng ###
def test_rng_is_concrete():
    assert_property_isconcrete(G_E_Phenotyping, "rng")

def test_rng_fget(ptprot):
    assert isinstance(ptprot.rng, (Generator,RandomState))

def test_rng_fset(ptprot):
    g = Generator(PCG64(12345))
    ptprot.rng = g
    assert id(ptprot.rng) == id(g)

def test_rng_fset_TypeError(ptprot):
    with pytest.raises(TypeError):
        ptprot.rng = object()

def test_rng_fdel(ptprot):
    with pytest.raises(AttributeError):
        del ptprot.rng

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_phenotype(ptprot, dpgmat, gpmod, nenv, nrep):
    # conduct true phenotyping
    df = ptprot.phenotype(dpgmat)

    # check if output is a pandas dataframe
    assert isinstance(df, pandas.DataFrame)

    # check that the number of rows == number of individuals
    assert len(df) == (dpgmat.ntaxa * nenv * nrep)

    # check that the number of columns == 4 + number of traits
    # four other columns are taxa, taxa_grp, env, rep
    assert len(df.columns) == (4 + gpmod.ntrait)

# TODO: finish these tests
def test_set_h2(ptprot, dpgmat):
    ptprot.set_h2(0.5, dpgmat)

def test_set_H2(ptprot, dpgmat):
    ptprot.set_H2(0.5, dpgmat)
