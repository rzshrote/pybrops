import numpy
import pytest

from pybrops.model.embvmat.DenseExpectedMaximumBreedingValueMatrix import DenseExpectedMaximumBreedingValueMatrix
from pybrops.model.embvmat.DenseExpectedMaximumBreedingValueMatrix import check_is_DenseExpectedMaximumBreedingValueMatrix
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.test.assert_python import assert_class_isconcrete, assert_method_isconcrete
from pybrops.test.assert_python import assert_module_public_api
from pybrops.test.assert_python import assert_classmethod_isconcrete
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

@pytest.fixture
def nprogeny():
    yield 20

@pytest.fixture
def nrep():
    yield 10

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
def pgmat_vrnt_xoprob(nmarker):
    out = numpy.random.uniform(0.0, 0.5, nmarker)
    yield out

@pytest.fixture
def pgmat(
        pgmat_mat, 
        pgmat_chrgrp, 
        pgmat_phypos, 
        pgmat_taxa, 
        pgmat_taxa_grp,
        pgmat_vrnt_xoprob,
    ):
    yield DensePhasedGenotypeMatrix(
        mat         = pgmat_mat,
        vrnt_chrgrp = pgmat_chrgrp,
        vrnt_phypos = pgmat_phypos,
        taxa        = pgmat_taxa,
        taxa_grp    = pgmat_taxa_grp,
        vrnt_xoprob = pgmat_vrnt_xoprob,
    )

############################################################
########## Expected Maximum Breeding Value Matrix ##########
############################################################
@pytest.fixture
def embvmat_mat(ntaxa, ntrait):
    out = numpy.random.random((ntaxa,ntrait))
    yield out

@pytest.fixture
def embvmat_location(ntrait):
    out = numpy.repeat(0.0, ntrait)
    yield out

@pytest.fixture
def embvmat_scale(ntrait):
    out = numpy.repeat(1.0, ntrait)
    yield out

@pytest.fixture
def embvmat_taxa(pgmat_taxa):
    yield pgmat_taxa

@pytest.fixture
def embvmat_taxa_grp(pgmat_taxa_grp):
    yield pgmat_taxa_grp

@pytest.fixture
def embvmat_trait(algmod_trait):
    yield algmod_trait

@pytest.fixture
def embvmat(
        embvmat_mat,
        embvmat_location,
        embvmat_scale,
        embvmat_taxa,
        embvmat_taxa_grp,
        embvmat_trait,
    ):
    out = DenseExpectedMaximumBreedingValueMatrix(
        mat      = embvmat_mat,
        location = embvmat_location,
        scale    = embvmat_scale,
        taxa     = embvmat_taxa,
        taxa_grp = embvmat_taxa_grp,
        trait    = embvmat_trait,
    )
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################
def test_DenseExpectedMaximumBreedingValueMatrix_module_documentation():
    import pybrops.model.embvmat.DenseExpectedMaximumBreedingValueMatrix
    assert_module_documentation(pybrops.model.embvmat.DenseExpectedMaximumBreedingValueMatrix)

def test_DenseExpectedMaximumBreedingValueMatrix_module_public_api():
    import pybrops.model.embvmat.DenseExpectedMaximumBreedingValueMatrix
    assert_module_public_api(pybrops.model.embvmat.DenseExpectedMaximumBreedingValueMatrix)

################################################################################
########################### Test class documentation ###########################
################################################################################
def test_DenseExpectedMaximumBreedingValueMatrix_is_concrete():
    assert_class_isconcrete(DenseExpectedMaximumBreedingValueMatrix)

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

### __init__
def test___init___is_concrete():
    assert_method_isconcrete(DenseExpectedMaximumBreedingValueMatrix, "__init__")

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

### from_gmod
def test_from_gmod_is_concrete():
    assert_classmethod_isconcrete(DenseExpectedMaximumBreedingValueMatrix, "from_gmod")

def test_from_gmod_TypeError(algmod, pgmat, nprogeny, nrep):
    with pytest.raises(TypeError):
        mat = DenseExpectedMaximumBreedingValueMatrix.from_gmod(None, pgmat, nprogeny, nrep)
    with pytest.raises(TypeError):
        mat = DenseExpectedMaximumBreedingValueMatrix.from_gmod(algmod, None, nprogeny, nrep)
    with pytest.raises(TypeError):
        mat = DenseExpectedMaximumBreedingValueMatrix.from_gmod(algmod, pgmat, None, nrep)
    with pytest.raises(TypeError):
        mat = DenseExpectedMaximumBreedingValueMatrix.from_gmod(algmod, pgmat, nprogeny, None)

def test_from_gmod(algmod, pgmat, nprogeny, nrep):
    # calculate EMBVs
    mat = DenseExpectedMaximumBreedingValueMatrix.from_gmod(algmod, pgmat, nprogeny, nrep)
    
    # assert output
    assert isinstance(mat, DenseExpectedMaximumBreedingValueMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################

### check_is_DenseExpectedMaximumBreedingValueMatrix
def test_check_is_DenseExpectedMaximumBreedingValueMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseExpectedMaximumBreedingValueMatrix)

def test_check_is_DenseExpectedMaximumBreedingValueMatrix(embvmat):
    with not_raises(TypeError):
        check_is_DenseExpectedMaximumBreedingValueMatrix(embvmat, "embvmat")
    with pytest.raises(TypeError):
        check_is_DenseExpectedMaximumBreedingValueMatrix(None, "embvmat")