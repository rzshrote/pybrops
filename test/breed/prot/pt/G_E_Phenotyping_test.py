from numbers import Integral
import os
from pathlib import Path
import numpy
from numpy.random import Generator, RandomState, PCG64
import pandas
import pytest
import copy
import h5py
from pybrops.test.assert_python import assert_classmethod_isconcrete, assert_property_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.breed.prot.pt.G_E_Phenotyping import G_E_Phenotyping
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

################################################################################
################################ Test fixtures #################################
################################################################################

@pytest.fixture
def nphase():
    yield 2

@pytest.fixture
def ntaxa():
    yield 20 # must be divisible by 4

@pytest.fixture
def ntaxa_grp():
    yield 4

@pytest.fixture
def nvrnt():
    yield 100 # must be divisible by 2

@pytest.fixture
def nchrom():
    yield 2

@pytest.fixture
def ntrait():
    yield 2

############################################################
######################## Genotypes #########################
############################################################
@pytest.fixture
def mat_int8(nphase, ntaxa, nvrnt):
    out = numpy.random.randint(0, 2, (nphase,ntaxa,nvrnt))
    out = out.astype("int8")
    yield out

@pytest.fixture
def mat_chrgrp(nvrnt, nchrom):
    out = numpy.repeat([1,2], nvrnt // nchrom)
    yield out

@pytest.fixture
def mat_phypos(nvrnt):
    out = numpy.arange(1, nvrnt+1)
    yield out

@pytest.fixture
def mat_genpos(nvrnt, nchrom):
    out = numpy.empty(nvrnt, dtype = float)
    nsnp = nvrnt // nchrom
    for i in range(nchrom):
        tmp = numpy.random.random(nsnp)
        tmp.sort()
        out[(nsnp*i):(nsnp*(i+1))] = tmp
    yield out

@pytest.fixture
def mat_taxa(ntaxa):
    out = numpy.array(["Line"+str(i).zfill(2) for i in range(ntaxa)], dtype = object)
    yield out

@pytest.fixture
def mat_taxa_grp(ntaxa, ntaxa_grp):
    out = numpy.repeat(list(range(ntaxa_grp)), ntaxa // ntaxa_grp)
    yield out

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
    out.group_taxa()
    out.group_vrnt()
    yield out

############################################################
###################### Genomic model #######################
############################################################
@pytest.fixture
def mat_beta(ntrait):
    out = numpy.random.random((1,ntrait))
    yield out

@pytest.fixture
def mat_u_misc():
    yield None

@pytest.fixture
def mat_u_a(nvrnt, ntrait):
    out = numpy.random.random((nvrnt,ntrait))
    yield out

@pytest.fixture
def trait(ntrait):
    out = numpy.array(["Trait"+str(i).zfill(2) for i in range(ntrait)], dtype = object)
    yield out

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
        var_err = var_err,
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(G_E_Phenotyping)

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

### __init__

def test___init___is_concrete():
    assert_method_isconcrete(G_E_Phenotyping, "__init__")

### __copy__

def test___copy___is_concrete():
    assert_method_isconcrete(G_E_Phenotyping, "__copy__")

### __deepcopy__

def test___deepcopy___is_concrete():
    assert_method_isconcrete(G_E_Phenotyping, "__deepcopy__")

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
############################# Test concrete methods ############################
################################################################################
        
### phenotype ###

def test_phenotype_is_concrete():
    assert_method_isconcrete(G_E_Phenotyping, "phenotype")

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

### set_h2 ###

def test_set_h2_is_concrete():
    assert_method_isconcrete(G_E_Phenotyping, "set_h2")

# TODO: finish this test
def test_set_h2(ptprot, dpgmat):
    ptprot.set_h2(0.5, dpgmat)

### set_H2 ###

def test_set_H2_is_concrete():
    assert_method_isconcrete(G_E_Phenotyping, "set_H2")

# TODO: finish this test
def test_set_H2(ptprot, dpgmat):
    ptprot.set_H2(0.5, dpgmat)

### to_hdf5

def test_to_hdf5_is_concrete():
    assert_method_isconcrete(G_E_Phenotyping, "to_hdf5")

def test_to_hdf5_str(ptprot):
    fp = "tmp.h5"
    with not_raises(Exception):
        ptprot.to_hdf5(fp)
    assert os.path.exists(fp)
    os.remove(fp)

def test_to_hdf5_Path(ptprot):
    fp = Path("tmp.h5")
    with not_raises(Exception):
        ptprot.to_hdf5(fp)
    assert os.path.exists(fp)
    os.remove(fp)

def test_to_hdf5_h5py_File(ptprot):
    fp = "tmp.h5"
    h5file = h5py.File(fp, "a")
    with not_raises(Exception):
        ptprot.to_hdf5(h5file)
    h5file.close()
    assert os.path.exists(fp)
    os.remove(fp)

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

### from_hdf5

def test_from_hdf5_is_concrete():
    assert_classmethod_isconcrete(G_E_Phenotyping, "from_hdf5")

def test_from_hdf5_str(ptprot, gpmod):
    fp = "tmp.h5"
    ptprot.to_hdf5(fp)
    out = G_E_Phenotyping.from_hdf5(fp, gpmod = gpmod)
    assert ptprot.nenv == out.nenv
    assert numpy.all(ptprot.nrep == out.nrep)
    assert numpy.all(ptprot.var_env == out.var_env)
    assert numpy.all(ptprot.var_rep == out.var_rep)
    assert numpy.all(ptprot.var_err == out.var_err)
    os.remove(fp)

def test_from_hdf5_Path(ptprot, gpmod):
    fp = Path("tmp.h5")
    ptprot.to_hdf5(fp)
    out = G_E_Phenotyping.from_hdf5(fp, gpmod = gpmod)
    assert ptprot.nenv == out.nenv
    assert numpy.all(ptprot.nrep == out.nrep)
    assert numpy.all(ptprot.var_env == out.var_env)
    assert numpy.all(ptprot.var_rep == out.var_rep)
    assert numpy.all(ptprot.var_err == out.var_err)
    os.remove(fp)

def test_from_hdf5_h5py_File(ptprot, gpmod):
    fp = Path("tmp.h5")
    ptprot.to_hdf5(fp)
    h5file = h5py.File(fp)
    out = G_E_Phenotyping.from_hdf5(h5file, gpmod = gpmod)
    assert ptprot.nenv == out.nenv
    assert numpy.all(ptprot.nrep == out.nrep)
    assert numpy.all(ptprot.var_env == out.var_env)
    assert numpy.all(ptprot.var_rep == out.var_rep)
    assert numpy.all(ptprot.var_err == out.var_err)
    h5file.close()
    os.remove(fp)

