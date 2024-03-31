import numpy
import pytest

from pybrops.model.wgebvmat.WeightedGenomicEstimatedBreedingValueMatrix import WeightedGenomicEstimatedBreedingValueMatrix
from pybrops.model.wgebvmat.WeightedGenomicEstimatedBreedingValueMatrix import check_is_WeightedGenomicEstimatedBreedingValueMatrix
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.test.assert_python import assert_class_isabstract, assert_module_public_api
from pybrops.test.assert_python import assert_classmethod_isabstract
from pybrops.test.assert_python import assert_function_isconcrete
from pybrops.test.assert_python import assert_module_documentation
from pybrops.test.assert_python import not_raises
from .common_fixtures import *

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
        beta        = algmod_beta,
        u_misc      = algmod_u_misc,
        u_a         = algmod_u_a,
        trait       = algmod_trait,
        model_name  = algmod_model_name,
        hyperparams = algmod_params,
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
        mat         = pgmat_mat,
        vrnt_chrgrp = pgmat_chrgrp,
        vrnt_phypos = pgmat_phypos,
        taxa        = pgmat_taxa,
        taxa_grp    = pgmat_taxa_grp,
    )

############################################################
########## Expected Maximum Breeding Value Matrix ##########
############################################################
@pytest.fixture
def wgebvmat():
    out = DummyWeightedGenomicEstimatedBreedingValueMatrix()
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################
def test_WeightedGenomicEstimatedBreedingValueMatrix_module_documentation():
    import pybrops.model.wgebvmat.WeightedGenomicEstimatedBreedingValueMatrix
    assert_module_documentation(pybrops.model.wgebvmat.WeightedGenomicEstimatedBreedingValueMatrix)

def test_WeightedGenomicEstimatedBreedingValueMatrix_module_public_api():
    import pybrops.model.wgebvmat.WeightedGenomicEstimatedBreedingValueMatrix
    assert_module_public_api(pybrops.model.wgebvmat.WeightedGenomicEstimatedBreedingValueMatrix)

################################################################################
########################### Test class documentation ###########################
################################################################################
def test_WeightedGenomicEstimatedBreedingValueMatrix_is_abstract():
    assert_class_isabstract(WeightedGenomicEstimatedBreedingValueMatrix)

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
########################## Test abstract classmethods ##########################
################################################################################

### from_algmod
def test_from_algmod_is_abstract():
    assert_classmethod_isabstract(WeightedGenomicEstimatedBreedingValueMatrix, "from_algmod")

################################################################################
######################### Test class utility functions #########################
################################################################################

### check_is_WeightedGenomicEstimatedBreedingValueMatrix
def test_check_is_WeightedGenomicEstimatedBreedingValueMatrix_is_concrete():
    assert_function_isconcrete(check_is_WeightedGenomicEstimatedBreedingValueMatrix)

def test_check_is_WeightedGenomicEstimatedBreedingValueMatrix(wgebvmat):
    with not_raises(TypeError):
        check_is_WeightedGenomicEstimatedBreedingValueMatrix(wgebvmat, "wgebvmat")
    with pytest.raises(TypeError):
        check_is_WeightedGenomicEstimatedBreedingValueMatrix(None, "wgebvmat")