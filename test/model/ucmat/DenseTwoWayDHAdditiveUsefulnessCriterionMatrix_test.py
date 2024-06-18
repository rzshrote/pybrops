from numbers import Integral
import numpy
import pytest

from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.model.ucmat.DenseTwoWayDHAdditiveUsefulnessCriterionMatrix import DenseTwoWayDHAdditiveUsefulnessCriterionMatrix
from pybrops.model.ucmat.DenseTwoWayDHAdditiveUsefulnessCriterionMatrix import check_is_DenseTwoWayDHAdditiveUsefulnessCriterionMatrix
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
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
####################### GEBV matrix ########################
############################################################

@pytest.fixture
def gebvmat(algmod, pgmat):
    out = algmod.gebv(pgmat)
    yield out

############################################################
####### Progeny Mean Estimated Breeding Value Matrix #######
############################################################

@pytest.fixture
def ucmat_mat(ntaxa, ntrait):
    out = numpy.random.random((ntaxa,ntaxa,ntrait))
    yield out

@pytest.fixture
def ucmat_taxa(ntaxa):
    out = numpy.array(["Line"+str(i).zfill(3) for i in range(ntaxa)], dtype=object)
    yield out

@pytest.fixture
def ucmat_taxa_grp(ntaxa, ngroup):
    out = numpy.random.randint(1, ngroup+1, ntaxa)
    out.sort()
    yield out

@pytest.fixture
def ucmat_trait(ntrait):
    out = numpy.array(["Trait"+str(i).zfill(3) for i in range(ntrait)], dtype=object)
    yield out

@pytest.fixture
def ucmat(
        ucmat_mat,
        ucmat_taxa,
        ucmat_taxa_grp,
        ucmat_trait,
    ):
    out = DenseTwoWayDHAdditiveUsefulnessCriterionMatrix(
        mat = ucmat_mat,
        taxa = ucmat_taxa,
        taxa_grp = ucmat_taxa_grp,
        trait = ucmat_trait,
    )
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix_module_documentation():
    import pybrops.model.ucmat.DenseTwoWayDHAdditiveUsefulnessCriterionMatrix
    assert_module_documentation(pybrops.model.ucmat.DenseTwoWayDHAdditiveUsefulnessCriterionMatrix)

def test_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix_module_public_api():
    import pybrops.model.ucmat.DenseTwoWayDHAdditiveUsefulnessCriterionMatrix
    assert_module_public_api(pybrops.model.ucmat.DenseTwoWayDHAdditiveUsefulnessCriterionMatrix)

################################################################################
############################ Test class attributes #############################
################################################################################

def test_DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix_is_concrete():
    assert_class_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix)

################################################################################
######################## Test abstract special methods #########################
################################################################################

################################################################################
######################## Test concrete special methods #########################
################################################################################

### __init__
def test___init___is_concrete():
    assert_method_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, "__init__")

### __copy__
def test___copy___is_concrete():
    assert_method_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, "__copy__")

### __deepcopy__
def test___deepcopy___is_concrete():
    assert_method_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, "__deepcopy__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
########################### Test concrete properties ###########################
################################################################################

### mat
def test_mat_is_concrete():
    assert_property_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, "mat")

def test_mat_fset_TypeError(ucmat):
    with pytest.raises(TypeError):
        ucmat.mat = object()

def test_mat_fset_ValueError(ucmat, ntaxa, ntrait):
    with pytest.raises(ValueError):
        ucmat.mat = numpy.random.random((ntaxa,))
    with pytest.raises(ValueError):
        ucmat.mat = numpy.random.random((ntaxa,ntrait))
    with pytest.raises(ValueError):
        ucmat.mat = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))

### square_taxa_axes
def test_square_taxa_axes_is_concrete():
    assert_property_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, "square_taxa_axes")

def test_square_taxa_axes_fget(ucmat):
    assert ucmat.square_taxa_axes == (0,1)

### trait_axis
def test_trait_axis_is_concrete():
    assert_property_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, "trait_axis")

def test_trait_axis_fget(ucmat):
    assert ucmat.trait_axis == 2

### epgc
def test_epgc_is_concrete():
    assert_property_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, "epgc")

def test_epgc_fget(ucmat):
    assert isinstance(ucmat.epgc, tuple)
    assert ucmat.epgc == (0.5, 0.5)

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################

### copy
def test_copy_is_concrete():
    assert_method_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, "copy")

### deepcopy
def test_deepcopy_is_concrete():
    assert_method_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, "deepcopy")

################################################################################
########################## Test abstract classmethods ##########################
################################################################################

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

### from_gmod
def test_from_gmod_is_concrete():
    assert_classmethod_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, "from_gmod")

def test_from_gmod(algmod, pgmat, gebvmat, ntaxa, ntrait):
    observed = DenseTwoWayDHAdditiveUsefulnessCriterionMatrix.from_gmod(
        gmod = algmod,
        pgmat = pgmat,
        nmating = 1,
        nprogeny = 80,
        nself = 0,
        gmapfn = HaldaneMapFunction(),
        upper_percentile = 0.05,
    )
    assert isinstance(observed, DenseTwoWayDHAdditiveUsefulnessCriterionMatrix)
    assert observed.mat_shape == (ntaxa,ntaxa,ntrait)
    # allocate expected matrix
    mean_mat = numpy.empty((ntaxa,ntaxa,ntrait), dtype = float)
    gebv = gebvmat.unscale()
    for i in range(ntaxa):
        for j in range(ntaxa):
            for k in range(ntrait):
                mean_mat[i,j,k] = 0.5 * (gebv[i,k] + gebv[j,k])
    assert numpy.all(observed.mat >= mean_mat)

### from_algmod
def test_from_algmod_is_concrete():
    assert_classmethod_isconcrete(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, "from_algmod")

def test_from_algmod(algmod, pgmat, gebvmat, ntaxa, ntrait):
    observed = DenseTwoWayDHAdditiveUsefulnessCriterionMatrix.from_algmod(
        algmod = algmod,
        pgmat = pgmat,
        nmating = 1,
        nprogeny = 80,
        nself = 0,
        gmapfn = HaldaneMapFunction(),
        upper_percentile = 0.05,
        mem = 1028,
    )
    assert isinstance(observed, DenseTwoWayDHAdditiveUsefulnessCriterionMatrix)
    assert observed.mat_shape == (ntaxa,ntaxa,ntrait)
    # allocate expected matrix
    mean_mat = numpy.empty((ntaxa,ntaxa,ntrait), dtype = float)
    gebv = gebvmat.unscale()
    for i in range(ntaxa):
        for j in range(ntaxa):
            for k in range(ntrait):
                mean_mat[i,j,k] = 0.5 * (gebv[i,k] + gebv[j,k])
    assert numpy.all(observed.mat >= mean_mat)

################################################################################
######################### Test abstract staticmethods ##########################
################################################################################

################################################################################
######################### Test concrete staticmethods ##########################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################

### check_is_DenseTwoWayDHAdditiveUsefulnessCriterionMatrix
def test_check_is_DenseTwoWayDHAdditiveUsefulnessCriterionMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseTwoWayDHAdditiveUsefulnessCriterionMatrix)

def test_check_is_DenseTwoWayDHAdditiveUsefulnessCriterionMatrix(ucmat):
    with not_raises(TypeError):
        check_is_DenseTwoWayDHAdditiveUsefulnessCriterionMatrix(ucmat, "ucmat")
    with pytest.raises(TypeError):
        check_is_DenseTwoWayDHAdditiveUsefulnessCriterionMatrix(None, "ucmat")