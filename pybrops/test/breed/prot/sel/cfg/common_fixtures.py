import random
import numpy
import pytest

from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

@pytest.fixture
def common_ncross():
    yield 10

@pytest.fixture
def common_nparent():
    yield 2

@pytest.fixture
def common_nmating():
    yield 1

@pytest.fixture
def common_nprogeny():
    yield 10

################################ Shape fixtures ################################
@pytest.fixture
def common_ntaxa():
    yield 200 # must be divisible by 4

@pytest.fixture
def common_nvrnt():
    yield 1000

@pytest.fixture
def common_ntrait():
    yield 2

@pytest.fixture
def common_nchr():
    yield 2

@pytest.fixture
def common_ngrp():
    yield 4

####################### Phased genotype matrix fixtures ########################

@pytest.fixture
def common_pgmat_mat(common_ntaxa, common_nvrnt):
    out = numpy.random.binomial(1, 0.5, (2,common_ntaxa,common_nvrnt)).astype('int8')
    yield out

@pytest.fixture
def common_taxa(common_ntaxa):
    out = numpy.array(["Taxa"+str(i).zfill(3) for i in range(1,common_ntaxa+1)], dtype=object)
    yield out

@pytest.fixture
def common_taxa_grp(common_ntaxa, common_ngrp):
    out = numpy.repeat(list(range(1,common_ngrp+1)), common_ntaxa // common_ngrp)
    yield out

@pytest.fixture
def common_vrnt_chrgrp(common_nvrnt, common_nchr):
    out = numpy.repeat(list(range(1,common_nchr+1)), common_nvrnt // common_nchr)
    yield out

@pytest.fixture
def common_vrnt_phypos(common_nvrnt):
    out = numpy.arange(1,common_nvrnt+1)
    yield out

@pytest.fixture
def common_vrnt_name(common_nvrnt):
    out = numpy.array(["SNP"+str(i).zfill(4) for i in range (1,common_nvrnt+1)], dtype=object)
    yield out

@pytest.fixture
def common_vrnt_genpos(common_nvrnt, common_nchr):
    l = []
    for i in range(common_nchr):
        tmp = numpy.random.random(common_nvrnt // common_nchr)
        tmp.sort()
        l.append(tmp)
    out = numpy.concatenate(l)
    yield out

@pytest.fixture
def common_vrnt_xoprob(common_nvrnt, common_nchr):
    l = []
    for i in range(common_nchr):
        tmp = numpy.random.uniform(0.0, 0.5, common_nvrnt // common_nchr)
        tmp[0] = 0.5
        l.append(tmp)
    out = numpy.concatenate(l)
    yield out

@pytest.fixture
def common_vrnt_hapgrp(common_nvrnt, common_nchr):
    l = []
    for i in range(common_nchr):
        tmp = numpy.random.randint(i*3+1, (i+1)*3+1, common_nvrnt // common_nchr)
        tmp.sort()
        l.append(tmp)
    out = numpy.concatenate(l)
    yield out

@pytest.fixture
def common_vrnt_hapref(common_nvrnt):
    out = numpy.array(random.choices("ACGT", k = common_nvrnt), dtype=object)
    yield out

@pytest.fixture
def common_vrnt_hapalt(common_vrnt_hapref):
    l = []
    for e in common_vrnt_hapref:
        if e == 'A':
            l.append(random.choice("CGT"))
        elif e == 'C':
            l.append(random.choice("AGT"))
        elif e == 'G':
            l.append(random.choice("ACT"))
        elif e == 'T':
            l.append(random.choice("ACG"))
    out = numpy.array(l, dtype=object)
    yield out

@pytest.fixture
def common_vrnt_mask(common_nvrnt):
    out = numpy.repeat(True, common_nvrnt)
    yield out

@pytest.fixture
def common_pgmat(
        common_pgmat_mat,
        common_taxa,
        common_taxa_grp,
        common_vrnt_chrgrp,
        common_vrnt_phypos,
        common_vrnt_name,
        common_vrnt_genpos,
        common_vrnt_xoprob, 
        common_vrnt_hapgrp, 
        common_vrnt_hapalt, 
        common_vrnt_hapref, 
        common_vrnt_mask
    ):
    out = DensePhasedGenotypeMatrix(
        mat = common_pgmat_mat,
        taxa = common_taxa,
        taxa_grp = common_taxa_grp,
        vrnt_chrgrp = common_vrnt_chrgrp,
        vrnt_phypos = common_vrnt_phypos,
        vrnt_name = common_vrnt_name,
        vrnt_genpos = common_vrnt_genpos,
        vrnt_xoprob = common_vrnt_xoprob, 
        vrnt_hapgrp = common_vrnt_hapgrp, 
        vrnt_hapalt = common_vrnt_hapalt, 
        vrnt_hapref = common_vrnt_hapref, 
        vrnt_mask = common_vrnt_mask 
    )
    out.group()
    yield out

##################### Cross configuration matrix fixtures ######################
@pytest.fixture
def common_xconfig(common_ncross, common_nparent, common_ntaxa):
    out = numpy.random.randint(0, common_ntaxa, (common_ncross, common_nparent))
    yield out
