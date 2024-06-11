import numpy
import pytest

from pybrops.breed.prot.sel.TwoWayParentalMeanGenomicEstimatedBreedingValueSelection import TwoWayParentalMeanGenomicEstimatedBreedingValueSubsetSelection
from pybrops.breed.prot.sel.prob.TwoWayParentalMeanGenomicEstimatedBreedingValueSelectionProblem import TwoWayParentalMeanGenomicEstimatedBreedingValueSubsetMateSelectionProblem
from pybrops.core.util.crossix import twowayix, twowayix_len
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.test.assert_python import assert_method_isconcrete

################################################################################
################################ Test fixtures #################################
################################################################################

###################### Shape fixtures ######################
@pytest.fixture
def ntaxa():
    yield 100

@pytest.fixture
def ntaxa_grp(ntaxa):
    out = ntaxa // 4
    yield out

@pytest.fixture
def ntrait():
    yield 3

@pytest.fixture
def nvrnt():
    yield 200

@pytest.fixture
def nvrnt_chrgrp(nvrnt):
    out = nvrnt // 20
    yield out

##################### Genotype Matrix ######################
@pytest.fixture
def gmat_mat(ntaxa,nvrnt):
    out = numpy.random.randint(0,3,(ntaxa,nvrnt)).astype("int8")
    yield out

@pytest.fixture
def gmat_taxa(ntaxa):
    iterator = ("Taxon"+str(i).zfill(3) for i in range(ntaxa))
    out = numpy.fromiter(iterator, dtype = object, count = ntaxa)
    yield out

@pytest.fixture
def gmat_taxa_grp(ntaxa,ntaxa_grp):
    out = numpy.repeat(numpy.arange(ntaxa_grp), ntaxa//ntaxa_grp)
    yield out

@pytest.fixture
def gmat_vrnt_chrgrp(nvrnt,nvrnt_chrgrp):
    out = numpy.repeat(numpy.arange(nvrnt_chrgrp), nvrnt//nvrnt_chrgrp)
    yield out

@pytest.fixture
def gmat_vrnt_phypos(nvrnt):
    out = numpy.random.randint(0,2**32,(nvrnt,))
    out.sort()
    yield out

@pytest.fixture
def gmat_vrnt_genpos(nvrnt):
    out = numpy.random.random((nvrnt,))
    out.sort()
    yield out

@pytest.fixture
def gmat(
        gmat_mat,
        gmat_taxa,
        gmat_taxa_grp,
        gmat_vrnt_chrgrp,
        gmat_vrnt_phypos,
        gmat_vrnt_genpos,
    ):
    out = DenseGenotypeMatrix(
        mat = gmat_mat,
        taxa = gmat_taxa,
        taxa_grp = gmat_taxa_grp,
        vrnt_chrgrp = gmat_vrnt_chrgrp,
        vrnt_phypos = gmat_vrnt_phypos,
        vrnt_genpos = gmat_vrnt_genpos,
    )
    yield out

############## Additive Linear Genomic Model ###############
@pytest.fixture
def algmod_beta(ntrait):
    out = numpy.random.random((1,ntrait))
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
    iterator = ("Trait"+str(i).zfill(2) for i in range(ntrait))
    out = numpy.fromiter(iterator, dtype = object, count = ntrait)
    yield out

@pytest.fixture
def algmod(
        algmod_beta,
        algmod_u_misc,
        algmod_u_a,
        algmod_trait,
    ):
    out = DenseAdditiveLinearGenomicModel(
        beta = algmod_beta,
        u_misc = algmod_u_misc,
        u_a = algmod_u_a,
        trait = algmod_trait,
    )
    yield out

####################### GEBV matrix ########################
@pytest.fixture
def gebvmat(algmod, gmat):
    out = algmod.gebv(gmat)
    yield out

########### 2-way pmGEBV mate selection problem ############
@pytest.fixture
def selprob_symab():
    yield True

@pytest.fixture
def selprob_mateab():
    yield "uniq"

@pytest.fixture
def selprob_xmap(ntaxa, selprob_symab, selprob_mateab):
    xmap = numpy.fromiter(
        iter = twowayix(ntaxa, selprob_symab, selprob_mateab),
        dtype = numpy.dtype((int, 2)),
        count = twowayix_len(ntaxa, selprob_symab, selprob_mateab)
    )
    yield xmap

@pytest.fixture
def selprob_unscale():
    yield True

# select 10 two-way crosses
@pytest.fixture
def selprob_ndecn():
    yield 10

@pytest.fixture
def selprob_decn_space(selprob_xmap):
    out = numpy.arange(len(selprob_xmap))
    yield out

@pytest.fixture
def selprob_decn_space_lower(selprob_ndecn):
    out = numpy.repeat(0, selprob_ndecn)
    yield out

@pytest.fixture
def selprob_decn_space_upper(selprob_ndecn,selprob_xmap):
    out = numpy.repeat(len(selprob_xmap)-1, selprob_ndecn)
    yield out

@pytest.fixture
def selprob(
        gmat,
        algmod,
        selprob_xmap,
        selprob_unscale,
        selprob_ndecn,
        selprob_decn_space,
        selprob_decn_space_lower,
        selprob_decn_space_upper,
        ntrait
    ):
    out = TwoWayParentalMeanGenomicEstimatedBreedingValueSubsetMateSelectionProblem.from_gmat_gpmod_xmap(
        gmat = gmat,
        gpmod = algmod,
        xmap = selprob_xmap,
        unscale = selprob_unscale,
        ndecn = selprob_ndecn,
        decn_space = selprob_decn_space,
        decn_space_lower = selprob_decn_space_lower,
        decn_space_upper = selprob_decn_space_upper,
        nobj = ntrait,
    )
    yield out

########### 2-way pmGEBV mate selection problem ############
@pytest.fixture
def selprot(
        ntrait,
    ):
    out = TwoWayParentalMeanGenomicEstimatedBreedingValueSubsetSelection(
        ntrait = ntrait,
        unscale = True,
        symab = True,
        mateab = "uniq",
        ncross = 20,
        nparent = 2,
        nmating = 1,
        nprogeny = 80,
        nobj = 3,
    )
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################

################################################################################
############################ Test class attributes #############################
################################################################################

################################################################################
######################## Test abstract special methods #########################
################################################################################

################################################################################
######################## Test concrete special methods #########################
################################################################################

### __init__
def test___init___is_concrete():
    assert_method_isconcrete(TwoWayParentalMeanGenomicEstimatedBreedingValueSubsetSelection, "__init__")

def test___init__(selprot):
    assert isinstance(selprot, TwoWayParentalMeanGenomicEstimatedBreedingValueSubsetSelection)

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################## Test abstract classmethods ##########################
################################################################################

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

################################################################################
######################### Test abstract staticmethods ##########################
################################################################################

################################################################################
######################### Test concrete staticmethods ##########################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
