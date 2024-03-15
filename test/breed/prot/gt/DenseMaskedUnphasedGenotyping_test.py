import numpy
import pytest

from pybrops.breed.prot.gt.DenseMaskedUnphasedGenotyping import DenseMaskedUnphasedGenotyping
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_class_documentation

from pybrops.breed.prot.gt.GenotypingProtocol import GenotypingProtocol

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
def mat(nphase, ntaxa, nvrnt):
    out = numpy.random.randint(0, 2, (nphase,ntaxa,nvrnt))
    out = out.astype("int8")
    yield out

@pytest.fixture
def taxa(ntaxa):
    out = numpy.array(["Line"+str(i).zfill(2) for i in range(ntaxa)], dtype = object)
    yield out

@pytest.fixture
def taxa_grp(ntaxa, ntaxa_grp):
    out = numpy.repeat(list(range(ntaxa_grp)), ntaxa // ntaxa_grp)
    yield out

@pytest.fixture
def vrnt_chrgrp(nvrnt, nchrom):
    out = numpy.repeat([1,2], nvrnt // nchrom)
    yield out

@pytest.fixture
def vrnt_phypos(nvrnt):
    out = numpy.arange(1, nvrnt+1)
    yield out

@pytest.fixture
def vrnt_genpos(nvrnt, nchrom):
    out = numpy.empty(nvrnt, dtype = float)
    nsnp = nvrnt // nchrom
    for i in range(nchrom):
        tmp = numpy.random.random(nsnp)
        tmp.sort()
        out[(nsnp*i):(nsnp*(i+1))] = tmp
    yield out

@pytest.fixture
def vrnt_mask(nvrnt):
    out = numpy.random.randint(0,2,nvrnt)
    out = out.astype(bool)
    yield out

@pytest.fixture
def nmask(vrnt_mask):
    out = numpy.sum(vrnt_mask)
    yield out

@pytest.fixture
def pgmat(mat, taxa, taxa_grp, vrnt_chrgrp, vrnt_phypos, vrnt_genpos, vrnt_mask):
    out = DensePhasedGenotypeMatrix(
        mat = mat,
        taxa = taxa,
        taxa_grp = taxa_grp,
        vrnt_chrgrp = vrnt_chrgrp,
        vrnt_phypos = vrnt_phypos,
        vrnt_genpos = vrnt_genpos,
        vrnt_mask = vrnt_mask,
    )
    out.group_taxa()
    out.group_vrnt()
    yield out

@pytest.fixture
def pgmat_nomask(mat, taxa, taxa_grp, vrnt_chrgrp, vrnt_phypos, vrnt_genpos):
    out = DensePhasedGenotypeMatrix(
        mat = mat,
        taxa = taxa,
        taxa_grp = taxa_grp,
        vrnt_chrgrp = vrnt_chrgrp,
        vrnt_phypos = vrnt_phypos,
        vrnt_genpos = vrnt_genpos,
    )
    out.group_taxa()
    out.group_vrnt()
    yield out

############################################################
################### Genotyping Protocol ####################
############################################################

@pytest.fixture
def gtprot():
    yield DenseMaskedUnphasedGenotyping()

################################################################################
############################ Test class documentation ##########################
################################################################################
def test_class_docstring():
    assert_class_documentation(DenseMaskedUnphasedGenotyping)

################################################################################
######################### Test concrete special methods ########################
################################################################################

### __init__

def test___init___is_concrete():
    assert_method_isconcrete(DenseMaskedUnphasedGenotyping, "__init__")

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################

### genotype

def test_genotype_is_concrete():
    assert_method_isconcrete(DenseMaskedUnphasedGenotyping, "genotype")

def test_genotype_TypeError(gtprot):
    with pytest.raises(TypeError):
        gtprot.genotype(object())

def test_genotype_mask(gtprot, pgmat, nmask):
    gmat = gtprot.genotype(pgmat)
    assert isinstance(gmat, GenotypeMatrix)
    assert gmat.ntaxa == pgmat.ntaxa
    assert gmat.nvrnt == nmask
    assert gmat.ploidy == pgmat.ploidy

def test_genotype_nomask(gtprot, pgmat_nomask, nvrnt):
    gmat = gtprot.genotype(pgmat_nomask)
    assert isinstance(gmat, GenotypeMatrix)
    assert gmat.ntaxa == pgmat_nomask.ntaxa
    assert gmat.nvrnt == nvrnt
    assert gmat.ploidy == pgmat_nomask.ploidy
