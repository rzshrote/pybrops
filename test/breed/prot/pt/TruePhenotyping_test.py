import os
from pathlib import Path
import numpy
import pandas
import pytest
import copy
import h5py
from pybrops.test.assert_python import assert_class_documentation, assert_classmethod_isconcrete, not_raises
from pybrops.test.assert_python import assert_method_isconcrete

from pybrops.breed.prot.pt.TruePhenotyping import TruePhenotyping
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
##################### TruePhenotyping ######################
############################################################
@pytest.fixture
def ptprot(gpmod):
    yield TruePhenotyping(
        gpmod = gpmod
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(TruePhenotyping)

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

### __init__

def test___init___is_concrete():
    assert_method_isconcrete(TruePhenotyping, "__init__")

### __copy__

def test___copy___is_concrete():
    assert_method_isconcrete(TruePhenotyping, "__copy__")

### __deepcopy__

def test___deepcopy___is_concrete():
    assert_method_isconcrete(TruePhenotyping, "__deepcopy__")

################################################################################
############################ Test Class Properties #############################
################################################################################

### gpmod

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

### var_err

def test_var_err_fget(ptprot):
    assert numpy.all(ptprot.var_err == 0.0)

def test_var_err_fset(ptprot):
    with pytest.raises(AttributeError):
        ptprot.var_err = 0.5

def test_var_err_fdel(ptprot):
    with pytest.raises(AttributeError):
        del ptprot.var_err

################################################################################
############################# Test concrete methods ############################
################################################################################

### phenotype

def test_phenotype_is_concrete():
    assert_method_isconcrete(TruePhenotyping, "phenotype")

def test_phenotype(ptprot, dpgmat, gpmod):
    # conduct true phenotyping
    df = ptprot.phenotype(dpgmat)

    # check if output is a pandas dataframe
    assert isinstance(df, pandas.DataFrame)

    # check that the number of rows == number of individuals
    assert len(df) == dpgmat.ntaxa

    # check that the number of columns == 2 + number of traits
    # two other columns are taxa and taxa_grp
    assert len(df.columns) == (2 + gpmod.ntrait)

### set_h2

def test_set_h2_is_concrete():
    assert_method_isconcrete(TruePhenotyping, "set_h2")

def test_set_h2(ptprot, dpgmat):
    with pytest.raises(AttributeError):
        ptprot.set_h2(0.5, dpgmat)

### set_H2

def test_set_H2_is_concrete():
    assert_method_isconcrete(TruePhenotyping, "set_H2")

def test_set_H2(ptprot, dpgmat):
    with pytest.raises(AttributeError):
        ptprot.set_H2(0.5, dpgmat)

### to_hdf5

def test_to_hdf5_is_concrete():
    assert_method_isconcrete(TruePhenotyping, "to_hdf5")

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
    assert_classmethod_isconcrete(TruePhenotyping, "from_hdf5")

def test_from_hdf5_str(ptprot, gpmod):
    fp = "tmp.h5"
    ptprot.to_hdf5(fp)
    out = TruePhenotyping.from_hdf5(fp, gpmod = gpmod)
    assert numpy.all(ptprot.var_err == out.var_err)
    os.remove(fp)

def test_from_hdf5_Path(ptprot, gpmod):
    fp = Path("tmp.h5")
    ptprot.to_hdf5(fp)
    out = TruePhenotyping.from_hdf5(fp, gpmod = gpmod)
    assert numpy.all(ptprot.var_err == out.var_err)
    os.remove(fp)

def test_from_hdf5_h5py_File(ptprot, gpmod):
    fp = Path("tmp.h5")
    ptprot.to_hdf5(fp)
    h5file = h5py.File(fp)
    out = TruePhenotyping.from_hdf5(h5file, gpmod = gpmod)
    assert numpy.all(ptprot.var_err == out.var_err)
    h5file.close()
    os.remove(fp)

