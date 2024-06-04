from numbers import Integral
import numpy
import pytest

from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.model.pmgebvmat.DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix import DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix
from pybrops.model.pmgebvmat.DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix import check_is_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
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

############################################################
############## Additive Linear Genomic Model ###############
############################################################

@pytest.fixture
def algmod_beta(nfixed, ntrait):
    out = numpy.random.random((nfixed,ntrait))
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
##################### Genotype Matrix ######################
############################################################

@pytest.fixture
def gmat_mat(ntaxa, nvrnt):
    out = numpy.random.binomial(2, 0.5, (ntaxa,nvrnt)).astype('int8')
    yield out

@pytest.fixture
def gmat_taxa(ntaxa):
    out = numpy.array(["Line"+str(i).zfill(3) for i in range(ntaxa)], dtype=object)
    yield out

@pytest.fixture
def gmat_taxa_grp(ntaxa, ngroup):
    out = numpy.random.randint(1, ngroup+1, ntaxa)
    out.sort()
    yield out

@pytest.fixture
def gmat(
        gmat_mat,
        gmat_taxa,
        gmat_taxa_grp,
    ):
    out = DenseGenotypeMatrix(
        mat = gmat_mat,
        taxa = gmat_taxa,
        taxa_grp = gmat_taxa_grp,
    )
    yield out

############################################################
####################### GEBV matrix ########################
############################################################

@pytest.fixture
def gebvmat(algmod, gmat):
    out = algmod.gebv(gmat)
    yield out

############################################################
####### Progeny Mean Estimated Breeding Value Matrix #######
############################################################

@pytest.fixture
def pmgebvmat_mat(ntaxa, ntrait):
    out = numpy.random.random((ntaxa,ntaxa,ntrait))
    yield out

@pytest.fixture
def pmgebvmat_location(ntrait):
    out = numpy.repeat(0.0, ntrait)
    yield out

@pytest.fixture
def pmgebvmat_scale(ntrait):
    out = numpy.repeat(1.0, ntrait)
    yield out

@pytest.fixture
def pmgebvmat_taxa(ntaxa):
    out = numpy.array(["Line"+str(i).zfill(3) for i in range(ntaxa)], dtype=object)
    yield out

@pytest.fixture
def pmgebvmat_taxa_grp(ntaxa, ngroup):
    out = numpy.random.randint(1, ngroup+1, ntaxa)
    out.sort()
    yield out

@pytest.fixture
def pmgebvmat_trait(ntrait):
    out = numpy.array(["Trait"+str(i).zfill(3) for i in range(ntrait)], dtype=object)
    yield out

@pytest.fixture
def pmgebvmat(
        pmgebvmat_mat,
        pmgebvmat_location,
        pmgebvmat_scale,
        pmgebvmat_taxa,
        pmgebvmat_taxa_grp,
        pmgebvmat_trait,
    ):
    out = DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix(
        mat = pmgebvmat_mat,
        location = pmgebvmat_location,
        scale = pmgebvmat_scale,
        taxa = pmgebvmat_taxa,
        taxa_grp = pmgebvmat_taxa_grp,
        trait = pmgebvmat_trait,
    )
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix_module_documentation():
    import pybrops.model.pmgebvmat.DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix
    assert_module_documentation(pybrops.model.pmgebvmat.DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix)

def test_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix_module_public_api():
    import pybrops.model.pmgebvmat.DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix
    assert_module_public_api(pybrops.model.pmgebvmat.DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix)

################################################################################
############################ Test class attributes #############################
################################################################################

def test_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix_is_concrete():
    assert_class_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix)

################################################################################
######################## Test abstract special methods #########################
################################################################################

################################################################################
######################## Test concrete special methods #########################
################################################################################

### __init__
def test___init___is_concrete():
    assert_method_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "__init__")

### __copy__
def test___copy___is_concrete():
    assert_method_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "__copy__")

### __deepcopy__
def test___deepcopy___is_concrete():
    assert_method_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "__deepcopy__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
########################### Test concrete properties ###########################
################################################################################

### mat
def test_mat_is_concrete():
    assert_property_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "mat")

def test_mat_fset_TypeError(pmgebvmat):
    with pytest.raises(TypeError):
        pmgebvmat.mat = object()

def test_mat_fset_ValueError(pmgebvmat, ntaxa, ntrait):
    with pytest.raises(ValueError):
        pmgebvmat.mat = numpy.random.random((ntaxa,))
    with pytest.raises(ValueError):
        pmgebvmat.mat = numpy.random.random((ntaxa,ntrait))
    with pytest.raises(ValueError):
        pmgebvmat.mat = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

### square_taxa_axes
def test_square_taxa_axes_is_concrete():
    assert_property_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "square_taxa_axes")

def test_square_taxa_axes_fget(pmgebvmat):
    assert pmgebvmat.square_taxa_axes == (0,1)

### trait_axis
def test_trait_axis_is_concrete():
    assert_property_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "trait_axis")

def test_trait_axis_fget(pmgebvmat):
    assert pmgebvmat.trait_axis == 2

### location
def test_location_is_concrete():
    assert_property_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "location")

def test_location_fset_TypeError(pmgebvmat):
    with pytest.raises(TypeError):
        pmgebvmat.location = object()

def test_location_fset_ValueError(pmgebvmat, ntrait):
    with pytest.raises(ValueError):
        pmgebvmat.location = numpy.random.random((ntrait+1,))
    with pytest.raises(ValueError):
        pmgebvmat.location = numpy.random.random((ntrait,ntrait))

### scale
def test_scale_is_concrete():
    assert_property_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "scale")

def test_scale_fset_TypeError(pmgebvmat):
    with pytest.raises(TypeError):
        pmgebvmat.scale = object()

def test_scale_fset_ValueError(pmgebvmat, ntrait):
    with pytest.raises(ValueError):
        pmgebvmat.scale = numpy.random.random((ntrait+1,))
    with pytest.raises(ValueError):
        pmgebvmat.scale = numpy.random.random((ntrait,ntrait))

### epgc
def test_epgc_is_concrete():
    assert_property_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "epgc")

def test_epgc_fget(pmgebvmat):
    assert isinstance(pmgebvmat.epgc, tuple)
    assert pmgebvmat.epgc == (0.5, 0.5)

### nfemale
def test_nfemale_is_concrete():
    assert_property_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "nfemale")

def test_nfemale(pmgebvmat, ntaxa):
    assert isinstance(pmgebvmat.nfemale, Integral)
    assert pmgebvmat.nfemale == ntaxa

### female_axis
def test_female_axis_is_concrete():
    assert_property_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "female_axis")

def test_female_axis(pmgebvmat):
    assert isinstance(pmgebvmat.female_axis, Integral)
    assert pmgebvmat.female_axis == 0

### nmale
def test_nmale_is_concrete():
    assert_property_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "nmale")

def test_nmale(pmgebvmat, ntaxa):
    assert isinstance(pmgebvmat.nmale, Integral)
    assert pmgebvmat.nmale == ntaxa

### male_axis
def test_male_axis_is_concrete():
    assert_property_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "male_axis")

def test_male_axis(pmgebvmat):
    assert isinstance(pmgebvmat.male_axis, Integral)
    assert pmgebvmat.male_axis == 1

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################

### copy
def test_copy_is_concrete():
    assert_method_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "copy")

### deepcopy
def test_deepcopy_is_concrete():
    assert_method_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "deepcopy")

################################################################################
########################## Test abstract classmethods ##########################
################################################################################

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

### from_gmod
def test_from_gmod_is_concrete():
    assert_classmethod_isconcrete(DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix, "from_gmod")

def test_from_gmod(algmod, gmat, gebvmat, ntaxa, ntrait):
    observed = DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix.from_gmod(
        gmod = algmod,
        gmat = gmat,
    )
    assert isinstance(observed, DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix)
    assert observed.mat_shape == (ntaxa,ntaxa,ntrait)
    # allocate expected matrix
    expected_mat = numpy.empty((ntaxa,ntaxa,ntrait), dtype = float)
    gebv = gebvmat.mat
    for i in range(ntaxa):
        for j in range(ntaxa):
            for k in range(ntrait):
                expected_mat[i,j,k] = 0.5 * (gebv[i,k] + gebv[j,k])
    assert numpy.allclose(observed.mat, expected_mat)

################################################################################
######################### Test abstract staticmethods ##########################
################################################################################

################################################################################
######################### Test concrete staticmethods ##########################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################

### check_is_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix
def test_check_is_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix)

def test_check_is_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix(pmgebvmat):
    with not_raises(TypeError):
        check_is_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix(pmgebvmat, "pmgebvmat")
    with pytest.raises(TypeError):
        check_is_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix(None, "pmgebvmat")