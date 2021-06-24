import pytest
import numpy

from pybropt.popgen.hmat import DenseHaplotypeMatrix
from pybropt.popgen.hmat import is_DenseHaplotypeMatrix
from pybropt.popgen.hmat import check_is_DenseHaplotypeMatrix
from pybropt.popgen.hmat import cond_check_is_DenseHaplotypeMatrix

################################################################################
################################## Genotypes ###################################
################################################################################
@pytest.fixture
def mat_int8():
    yield numpy.int8(
        [[1, 0, 1, 0, 2, 2, 0, 1, 0, 1, 1, 2, 0, 1, 1, 1],
         [2, 2, 0, 1, 1, 2, 0, 0, 0, 2, 0, 0, 0, 1, 1, 1],
         [1, 1, 1, 0, 1, 1, 0, 1, 2, 1, 2, 2, 0, 1, 2, 1],
         [1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 0, 0, 0, 1, 1],
         [1, 1, 2, 1, 1, 0, 1, 0, 1, 1, 1, 2, 1, 0, 1, 1],
         [1, 1, 2, 1, 1, 0, 1, 1, 0, 2, 1, 1, 1, 2, 0, 2],
         [1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 2, 0, 1, 2],
         [1, 0, 0, 1, 2, 1, 1, 1, 2, 0, 1, 2, 2, 1, 1, 1],
         [1, 0, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 2, 2, 1, 1],
         [0, 1, 2, 1, 2, 1, 1, 0, 0, 2, 0, 1, 0, 1, 1, 2],
         [1, 1, 1, 0, 2, 0, 2, 1, 1, 1, 2, 1, 2, 0, 0, 0],
         [0, 2, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 2, 2, 1],
         [2, 1, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 0, 1],
         [2, 2, 1, 2, 0, 0, 1, 0, 1, 2, 1, 1, 0, 2, 1, 1],
         [2, 0, 2, 2, 0, 0, 0, 0, 1, 2, 1, 1, 1, 1, 1, 2],
         [2, 0, 1, 0, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2, 0, 1],
         [0, 0, 0, 0, 0, 2, 2, 1, 1, 1, 1, 1, 0, 1, 0, 1],
         [1, 2, 1, 1, 1, 0, 1, 1, 2, 0, 0, 1, 1, 1, 0, 1],
         [1, 1, 1, 1, 2, 1, 1, 0, 1, 0, 2, 1, 1, 1, 2, 1],
         [0, 1, 0, 1, 1, 0, 0, 1, 2, 1, 1, 2, 0, 0, 2, 0],
         [2, 1, 2, 2, 1, 1, 0, 1, 1, 1, 2, 2, 1, 0, 1, 2],
         [1, 1, 0, 0, 1, 0, 1, 0, 1, 2, 0, 1, 1, 1, 2, 1],
         [0, 1, 1, 2, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0],
         [1, 2, 0, 2, 1, 1, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1],
         [1, 1, 1, 0, 0, 1, 0, 2, 1, 0, 1, 1, 1, 2, 1, 2],
         [0, 1, 1, 1, 0, 1, 1, 2, 0, 1, 0, 1, 2, 0, 1, 1],
         [0, 2, 1, 1, 0, 1, 1, 1, 0, 1, 1, 2, 0, 1, 1, 1],
         [1, 2, 2, 1, 1, 0, 2, 1, 1, 1, 0, 1, 1, 2, 2, 1],
         [1, 2, 1, 0, 1, 0, 1, 1, 2, 0, 2, 1, 2, 0, 1, 1],
         [1, 0, 0, 2, 0, 1, 0, 0, 1, 1, 2, 1, 1, 1, 1, 0]]
    )

@pytest.fixture
def mat_taxa():
    yield numpy.object_(["Line"+str(i).zfill(2) for i in range(30)])

@pytest.fixture
def mat_taxa_grp():
    yield numpy.repeat([0,1,2], 10)

@pytest.fixture
def mat_chrgrp():
    yield numpy.int64([1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2])

@pytest.fixture
def mat_phypos():
    yield numpy.arange(16)

@pytest.fixture
def mat_name():
    yield numpy.object_([chr(i) for i in range(ord('a'), ord('a')+16)])

@pytest.fixture
def mat_genpos():
    yield numpy.linspace(0,1,16)

@pytest.fixture
def mat_xoprob():
    yield numpy.repeat(0.1, 16)

@pytest.fixture
def mat_hapgrp():
    yield numpy.repeat([i for i in range(4)], 4)

@pytest.fixture
def mat_hapalt():
    yield numpy.object_(['A','T','C','G']*4)

@pytest.fixture
def mat_hapref():
    yield numpy.object_(['T','A','G','C']*4)

@pytest.fixture
def mat_mask():
    yield numpy.repeat([True, False], 8)

@pytest.fixture
def mat_ploidy():
    yield 2

@pytest.fixture
def dhmat(mat_int8, mat_taxa, mat_taxa_grp, mat_chrgrp, mat_phypos, mat_name, mat_genpos, mat_xoprob, mat_hapgrp, mat_hapalt, mat_hapref, mat_mask, mat_ploidy):
    yield DenseHaplotypeMatrix(
        mat = mat_int8,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp,
        vrnt_chrgrp = mat_chrgrp,
        vrnt_phypos = mat_phypos,
        vrnt_name = mat_name,
        vrnt_genpos = mat_genpos,
        vrnt_xoprob = mat_xoprob,
        vrnt_hapgrp = mat_hapgrp,
        vrnt_hapalt = mat_hapalt,
        vrnt_hapref = mat_hapref,
        vrnt_mask = mat_mask,
        ploidy = mat_ploidy
    )

################################################################################
#################################### Tests #####################################
################################################################################

############################################################
###################### Function Tests ######################
############################################################
def test_init(dhmat):
    assert dhmat is not None

def test_is_DenseHaplotypeMatrix(dhmat):
    assert is_DenseHaplotypeMatrix(dhmat)

def test_check_is_DenseHaplotypeMatrix(dhmat):
    check_is_DenseHaplotypeMatrix(dhmat, "dhmat")
    with pytest.raises(TypeError):
        check_is_DenseHaplotypeMatrix(None, "None")

def test_cond_check_is_DenseHaplotypeMatrix(dhmat):
    cond_check_is_DenseHaplotypeMatrix(dhmat, "dhmat")
    cond_check_is_DenseHaplotypeMatrix(None, "None")
    with pytest.raises(TypeError):
        cond_check_is_DenseHaplotypeMatrix("incorrect type", "str")

############################################################
###################### Property Tests ######################
############################################################

########################################
####### DenseHaplotypeMatrix.mat #######
########################################
def test_mat_fget(dhmat, mat_int8):
    assert numpy.all(dhmat.mat == mat_int8)

def test_mat_fset(dhmat, mat_int8):
    with pytest.raises(TypeError):
        dhmat.mat = "string type"
    with pytest.raises(TypeError):
        dhmat.mat = numpy.int16(mat_int8)
    with pytest.raises(ValueError):
        dhmat.mat = mat_int8[0]

def test_mat_fdel(dhmat):
    del dhmat.mat
    with pytest.raises(AttributeError):
        dhmat.mat

########################################
##### DenseHaplotypeMatrix.ploidy ######
########################################
def test_ploidy_fget(dhmat, mat_ploidy):
    assert dhmat.ploidy == mat_ploidy

def test_ploidy_fset(dhmat, mat_ploidy):
    with pytest.raises(AttributeError):
        dhmat.ploidy = mat_ploidy

def test_ploidy_fdel(dhmat, mat_ploidy):
    with pytest.raises(AttributeError):
        del dhmat.ploidy

########################################
##### DenseHaplotypeMatrix.nphase ######
########################################
def test_nphase_fget(dhmat):
    assert dhmat.nphase == 0

def test_nphase_fset(dhmat):
    with pytest.raises(AttributeError):
        dhmat.nphase = 0

def test_nphase_fdel(dhmat):
    with pytest.raises(AttributeError):
        del dhmat.nphase

########################################
####### DenseHaplotypeMatrix.mat_format #######
########################################
def test_mat_format_fget(dhmat):
    assert isinstance(dhmat.mat_format, str)

########################################
#### DenseHaplotypeMatrix.mat_ndim #####
########################################
def test_mat_ndim_fget(dhmat):
    assert dhmat.mat_ndim == 2

def test_mat_ndim_fset(dhmat):
    with pytest.raises(AttributeError):
        dhmat.mat_ndim = 2

def test_mat_ndim_fdel(dhmat):
    with pytest.raises(AttributeError):
        del dhmat.mat_ndim

########################################
####### DenseHaplotypeMatrix.taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.taxa_grp #######
########################################

########################################
####### DenseHaplotypeMatrix.ntaxa #######
########################################

########################################
####### DenseHaplotypeMatrix.taxa_axis #######
########################################

########################################
####### DenseHaplotypeMatrix.taxa_grp_name #######
########################################

########################################
####### DenseHaplotypeMatrix.taxa_grp_stix #######
########################################

########################################
####### DenseHaplotypeMatrix.taxa_grp_spix #######
########################################

########################################
####### DenseHaplotypeMatrix.taxa_grp_len #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_chrgrp #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_phypos #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_name #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_genpos #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_xoprob #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_hapgrp #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_hapalt #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_hapref #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_mask #######
########################################

########################################
####### DenseHaplotypeMatrix.nvrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_axis #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_chrgrp_name #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_chrgrp_stix #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_chrgrp_spix #######
########################################

########################################
####### DenseHaplotypeMatrix.vrnt_chrgrp_len #######
########################################

############################################################
####################### Method Tests #######################
############################################################

########################################
####### DenseHaplotypeMatrix.adjoin #######
########################################
def test_adjoin(dhmat, mat_int8, mat_taxa, mat_taxa_grp, mat_chrgrp, mat_phypos, mat_name, mat_genpos, mat_xoprob, mat_hapgrp, mat_hapalt, mat_hapref, mat_mask):
    a = dhmat.adjoin(dhmat, axis = 0)
    assert numpy.all(a.mat == numpy.concatenate([mat_int8, mat_int8], axis = 0))
    assert numpy.all(a.taxa == numpy.concatenate([mat_taxa, mat_taxa], axis = 0))
    assert numpy.all(a.taxa_grp == numpy.concatenate([mat_taxa_grp, mat_taxa_grp], axis = 0))
    assert numpy.all(a.vrnt_chrgrp == mat_chrgrp)
    assert numpy.all(a.vrnt_phypos == mat_phypos)
    assert numpy.all(a.vrnt_name == mat_name)
    assert numpy.all(a.vrnt_genpos == mat_genpos)
    assert numpy.all(a.vrnt_xoprob == mat_xoprob)
    assert numpy.all(a.vrnt_hapgrp == mat_hapgrp)
    assert numpy.all(a.vrnt_hapalt == mat_hapalt)
    assert numpy.all(a.vrnt_hapref == mat_hapref)
    assert numpy.all(a.vrnt_mask == mat_mask)
    a = dhmat.adjoin(dhmat, axis = 1)
    assert numpy.all(a.mat == numpy.concatenate([mat_int8, mat_int8], axis = 1))
    assert numpy.all(a.taxa == mat_taxa)
    assert numpy.all(a.taxa_grp == mat_taxa_grp)
    assert numpy.all(a.vrnt_chrgrp == numpy.concatenate([mat_chrgrp, mat_chrgrp], axis = 0))
    assert numpy.all(a.vrnt_phypos == numpy.concatenate([mat_phypos, mat_phypos], axis = 0))
    assert numpy.all(a.vrnt_name == numpy.concatenate([mat_name, mat_name], axis = 0))
    assert numpy.all(a.vrnt_genpos == numpy.concatenate([mat_genpos, mat_genpos], axis = 0))
    assert numpy.all(a.vrnt_xoprob == numpy.concatenate([mat_xoprob, mat_xoprob], axis = 0))
    assert numpy.all(a.vrnt_hapgrp == numpy.concatenate([mat_hapgrp, mat_hapgrp], axis = 0))
    assert numpy.all(a.vrnt_hapalt == numpy.concatenate([mat_hapalt, mat_hapalt], axis = 0))
    assert numpy.all(a.vrnt_hapref == numpy.concatenate([mat_hapref, mat_hapref], axis = 0))
    assert numpy.all(a.vrnt_mask == numpy.concatenate([mat_mask, mat_mask], axis = 0))


########################################
####### DenseHaplotypeMatrix.adjoin_taxa #######
########################################
def test_adjoin_taxa(dhmat, mat_int8, mat_taxa, mat_taxa_grp, mat_chrgrp, mat_phypos, mat_name, mat_genpos, mat_xoprob, mat_hapgrp, mat_hapalt, mat_hapref, mat_mask):
    a = dhmat.adjoin_taxa(dhmat)
    assert numpy.all(a.mat == numpy.concatenate([mat_int8, mat_int8], axis = 0))
    assert numpy.all(a.taxa == numpy.concatenate([mat_taxa, mat_taxa], axis = 0))
    assert numpy.all(a.taxa_grp == numpy.concatenate([mat_taxa_grp, mat_taxa_grp], axis = 0))
    assert numpy.all(a.vrnt_chrgrp == mat_chrgrp)
    assert numpy.all(a.vrnt_phypos == mat_phypos)
    assert numpy.all(a.vrnt_name == mat_name)
    assert numpy.all(a.vrnt_genpos == mat_genpos)
    assert numpy.all(a.vrnt_xoprob == mat_xoprob)
    assert numpy.all(a.vrnt_hapgrp == mat_hapgrp)
    assert numpy.all(a.vrnt_hapalt == mat_hapalt)
    assert numpy.all(a.vrnt_hapref == mat_hapref)
    assert numpy.all(a.vrnt_mask == mat_mask)

########################################
####### DenseHaplotypeMatrix.adjoin_vrnt #######
########################################
def test_adjoin_vrnt(dhmat, mat_int8, mat_taxa, mat_taxa_grp, mat_chrgrp, mat_phypos, mat_name, mat_genpos, mat_xoprob, mat_hapgrp, mat_hapalt, mat_hapref, mat_mask):
    a = dhmat.adjoin_vrnt(dhmat)
    assert numpy.all(a.mat == numpy.concatenate([mat_int8, mat_int8], axis = 1))
    assert numpy.all(a.taxa == mat_taxa)
    assert numpy.all(a.taxa_grp == mat_taxa_grp)
    assert numpy.all(a.vrnt_chrgrp == numpy.concatenate([mat_chrgrp, mat_chrgrp], axis = 0))
    assert numpy.all(a.vrnt_phypos == numpy.concatenate([mat_phypos, mat_phypos], axis = 0))
    assert numpy.all(a.vrnt_name == numpy.concatenate([mat_name, mat_name], axis = 0))
    assert numpy.all(a.vrnt_genpos == numpy.concatenate([mat_genpos, mat_genpos], axis = 0))
    assert numpy.all(a.vrnt_xoprob == numpy.concatenate([mat_xoprob, mat_xoprob], axis = 0))
    assert numpy.all(a.vrnt_hapgrp == numpy.concatenate([mat_hapgrp, mat_hapgrp], axis = 0))
    assert numpy.all(a.vrnt_hapalt == numpy.concatenate([mat_hapalt, mat_hapalt], axis = 0))
    assert numpy.all(a.vrnt_hapref == numpy.concatenate([mat_hapref, mat_hapref], axis = 0))
    assert numpy.all(a.vrnt_mask == numpy.concatenate([mat_mask, mat_mask], axis = 0))

########################################
####### DenseHaplotypeMatrix.delete #######
########################################

########################################
####### DenseHaplotypeMatrix.delete_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.delete_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.insert #######
########################################

########################################
####### DenseHaplotypeMatrix.insert_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.insert_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.select #######
########################################

########################################
####### DenseHaplotypeMatrix.select_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.select_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.concat #######
########################################

########################################
####### DenseHaplotypeMatrix.concat_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.concat_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.append #######
########################################

########################################
####### DenseHaplotypeMatrix.append_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.append_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.remove #######
########################################

########################################
####### DenseHaplotypeMatrix.remove_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.remove_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.incorp #######
########################################

########################################
####### DenseHaplotypeMatrix.incorp_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.incorp_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.lexsort #######
########################################

########################################
####### DenseHaplotypeMatrix.lexsort_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.lexsort_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.reorder #######
########################################

########################################
####### DenseHaplotypeMatrix.reorder_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.reorder_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.sort #######
########################################

########################################
####### DenseHaplotypeMatrix.sort_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.sort_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.group #######
########################################

########################################
####### DenseHaplotypeMatrix.group_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.group_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.is_grouped #######
########################################

########################################
####### DenseHaplotypeMatrix.is_grouped_taxa #######
########################################

########################################
####### DenseHaplotypeMatrix.is_grouped_vrnt #######
########################################

########################################
####### DenseHaplotypeMatrix.thcount #######
########################################

########################################
####### DenseHaplotypeMatrix.thfreq #######
########################################

########################################
####### DenseHaplotypeMatrix.hcount #######
########################################

########################################
####### DenseHaplotypeMatrix.hfreq #######
########################################

########################################
####### DenseHaplotypeMatrix.mhf #######
########################################

########################################
####### DenseHaplotypeMatrix.mehe #######
########################################

########################################
####### DenseHaplotypeMatrix.gtcount #######
########################################

########################################
####### DenseHaplotypeMatrix.gtfreq #######
########################################
