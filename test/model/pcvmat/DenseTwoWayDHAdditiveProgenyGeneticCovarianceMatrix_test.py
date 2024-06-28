from numbers import Integral
import numpy
import pytest

from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.model.pcvmat.DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix import DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, check_is_DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.test.assert_python import assert_class_isconcrete
from pybrops.test.assert_python import assert_classmethod_isconcrete
from pybrops.test.assert_python import assert_property_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_module_documentation
from pybrops.test.assert_python import assert_module_public_api
from pybrops.test.assert_python import not_raises

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
################# General shape parameters #################
############################################################

@pytest.fixture
def ngroup():
    yield 5

@pytest.fixture
def ntaxa(ngroup):
    yield ngroup * 20

@pytest.fixture
def ntrait():
    yield 3

@pytest.fixture
def nfixed():
    yield 1

@pytest.fixture
def nvrnt():
    yield 100

@pytest.fixture
def nphase():
    yield 2

@pytest.fixture
def nchrom():
    yield 5

############################################################
############## Additive Linear Genomic Model ###############
############################################################

@pytest.fixture
def algmod_mean(ntrait):
    out = numpy.repeat(0.0, ntrait)
    yield out

@pytest.fixture
def algmod_cov(ntrait):
    out = numpy.array([
        [ 1.0, -0.4,  0.8],
        [-0.4,  1.0, -0.6],
        [ 0.8, -0.6,  1.0],
    ])
    yield out

@pytest.fixture
def algmod_beta(nfixed, algmod_mean, algmod_cov):
    out = numpy.random.multivariate_normal(
        mean = algmod_mean,
        cov = algmod_cov,
        size = (nfixed,)
    )
    yield out

@pytest.fixture
def algmod_u_misc():
    yield None

@pytest.fixture
def algmod_u_a(nvrnt,ntrait):
    out = numpy.random.random((nvrnt,ntrait))
    yield out

@pytest.fixture
def algmod_trait(ntrait):
    out = numpy.array(["Trait"+str(i) for i in range(ntrait)], dtype = object)
    yield out

@pytest.fixture
def algmod_model_name():
    yield "dummy"

@pytest.fixture
def algmod_hyperparams():
    yield None

@pytest.fixture
def algmod(
        algmod_beta,
        algmod_u_misc,
        algmod_u_a,
        algmod_trait,
        algmod_model_name,
        algmod_hyperparams,
    ):
    out = DenseAdditiveLinearGenomicModel(
        beta = algmod_beta,
        u_misc = algmod_u_misc,
        u_a = algmod_u_a,
        trait = algmod_trait,
        model_name = algmod_model_name,
        hyperparams = algmod_hyperparams,
    )
    yield out

############################################################
################## Phased Genotype Matrix ##################
############################################################

@pytest.fixture
def pgmat_mat(nphase, ntaxa, nvrnt):
    out = numpy.random.binomial(1, 0.5, (nphase,ntaxa,nvrnt)).astype('int8')
    yield out

@pytest.fixture
def pgmat_taxa(ntaxa):
    out = numpy.array(["Line"+str(i).zfill(3) for i in range(ntaxa)], dtype=object)
    yield out

@pytest.fixture
def pgmat_taxa_grp(ntaxa, ngroup):
    out = numpy.random.randint(1, ngroup+1, ntaxa)
    out.sort()
    yield out

@pytest.fixture
def pgmat_vrnt_chrgrp(nvrnt, nchrom):
    out = numpy.repeat(numpy.arange(nvrnt//nchrom) + 1, nchrom)
    yield out

@pytest.fixture
def pgmat_vrnt_phypos(nvrnt):
    out = numpy.random.randint(0,2**32,nvrnt)
    out.sort()
    yield out

@pytest.fixture
def pgmat_vrnt_genpos(nvrnt):
    out = numpy.random.random(nvrnt)
    out.sort()
    yield out

@pytest.fixture
def pgmat(
        pgmat_mat,
        pgmat_taxa,
        pgmat_taxa_grp,
        pgmat_vrnt_chrgrp,
        pgmat_vrnt_phypos,
        pgmat_vrnt_genpos,
    ):
    out = DensePhasedGenotypeMatrix(
        mat = pgmat_mat,
        taxa = pgmat_taxa,
        taxa_grp = pgmat_taxa_grp,
        vrnt_chrgrp = pgmat_vrnt_chrgrp,
        vrnt_phypos = pgmat_vrnt_phypos,
        vrnt_genpos = pgmat_vrnt_genpos,
    )
    out.group_vrnt()
    yield out

############################################################
####### Progeny Mean Estimated Breeding Value Matrix #######
############################################################

@pytest.fixture
def pcvmat_mat(ntaxa, ntrait):
    out = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))
    yield out

@pytest.fixture
def pcvmat_taxa(ntaxa):
    out = numpy.array(["Line"+str(i).zfill(3) for i in range(ntaxa)], dtype=object)
    yield out

@pytest.fixture
def pcvmat_taxa_grp(ntaxa, ngroup):
    out = numpy.random.randint(1, ngroup+1, ntaxa)
    out.sort()
    yield out

@pytest.fixture
def pcvmat_trait(ntrait):
    out = numpy.array(["Trait"+str(i).zfill(3) for i in range(ntrait)], dtype=object)
    yield out

@pytest.fixture
def pcvmat(
        pcvmat_mat,
        pcvmat_taxa,
        pcvmat_taxa_grp,
        pcvmat_trait,
    ):
    out = DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix(
        mat = pcvmat_mat,
        taxa = pcvmat_taxa,
        taxa_grp = pcvmat_taxa_grp,
        trait = pcvmat_trait,
    )
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix_module_documentation():
    import pybrops.model.pcvmat.DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix
    assert_module_documentation(pybrops.model.pcvmat.DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix)

def test_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix_module_public_api():
    import pybrops.model.pcvmat.DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix
    assert_module_public_api(pybrops.model.pcvmat.DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix)

################################################################################
############################ Test class attributes #############################
################################################################################

def test_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix_is_concrete():
    assert_class_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix)

################################################################################
######################## Test abstract special methods #########################
################################################################################

################################################################################
######################## Test concrete special methods #########################
################################################################################

### __init__
def test___init___is_concrete():
    assert_method_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, "__init__")

### __copy__
def test___copy___is_concrete():
    assert_method_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, "__copy__")

### __deepcopy__
def test___deepcopy___is_concrete():
    assert_method_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, "__deepcopy__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
########################### Test concrete properties ###########################
################################################################################

### mat
def test_mat_is_concrete():
    assert_property_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, "mat")

def test_mat_fset_TypeError(pcvmat):
    with pytest.raises(TypeError):
        pcvmat.mat = object()

def test_mat_fset_ValueError(pcvmat, ntaxa, ntrait):
    with pytest.raises(ValueError):
        pcvmat.mat = numpy.random.random((ntaxa,))
    with pytest.raises(ValueError):
        pcvmat.mat = numpy.random.random((ntaxa,ntrait))
    with pytest.raises(ValueError):
        pcvmat.mat = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait,ntrait))

### square_taxa_axes
def test_square_taxa_axes_is_concrete():
    assert_property_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, "square_taxa_axes")

def test_square_taxa_axes_fget(pcvmat):
    assert pcvmat.square_taxa_axes == (0,1)

### trait_axis
def test_trait_axis_is_concrete():
    assert_property_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, "trait_axis")

def test_trait_axis_fget(pcvmat):
    assert pcvmat.trait_axis == 2

### epgc
def test_epgc_is_concrete():
    assert_property_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, "epgc")

def test_epgc_fget(pcvmat):
    assert isinstance(pcvmat.epgc, tuple)
    assert pcvmat.epgc == (0.5, 0.5)

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################

### copy
def test_copy_is_concrete():
    assert_method_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, "copy")

### deepcopy
def test_deepcopy_is_concrete():
    assert_method_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, "deepcopy")

################################################################################
########################## Test abstract classmethods ##########################
################################################################################

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

### from_gmod
def test_from_gmod_is_concrete():
    assert_classmethod_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, "from_gmod")

def test_from_gmod(algmod, pgmat, ntaxa, ntrait):
    observed = DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix.from_gmod(
        gmod = algmod,
        pgmat = pgmat,
        nmating = 1,
        nprogeny = 80,
        nself = 0,
        gmapfn = HaldaneMapFunction(),
    )
    assert isinstance(observed, DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix)
    assert observed.mat_shape == (ntaxa,ntaxa,ntrait,ntrait)

### from_algmod
def test_from_algmod_is_concrete():
    assert_classmethod_isconcrete(DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix, "from_algmod")

def test_from_algmod(algmod, pgmat, ntaxa, ntrait):
    observed = DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix.from_algmod(
        algmod = algmod,
        pgmat = pgmat,
        nmating = 1,
        nprogeny = 80,
        nself = 0,
        gmapfn = HaldaneMapFunction(),
        mem = 1028,
    )
    assert isinstance(observed, DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix)
    assert observed.mat_shape == (ntaxa,ntaxa,ntrait,ntrait)
    print(observed.mat[0,1,:,:])
    assert False

################################################################################
######################### Test abstract staticmethods ##########################
################################################################################

################################################################################
######################### Test concrete staticmethods ##########################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################

### check_is_DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix
def test_check_is_DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix)

def test_check_is_DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix(pcvmat):
    with not_raises(TypeError):
        check_is_DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix(pcvmat, "pcvmat")
    with pytest.raises(TypeError):
        check_is_DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix(None, "pcvmat")