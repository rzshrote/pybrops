import random
import numpy
import pytest

from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.core.random.prng import global_prng

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
def common_nconfig(common_ntaxa):
    yield (common_ntaxa * (common_ntaxa - 1)) // 2

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
def common_decn_space_xmap(common_nconfig, common_nparent, common_ntaxa):
    out = numpy.random.randint(0, common_ntaxa, (common_nconfig, common_nparent))
    yield out

### constructor parameters
@pytest.fixture
def common_ndecn_mate_subset(common_nconfig):
    out = common_nconfig // 2
    yield out

@pytest.fixture
def common_ndecn_mate_binary(common_nconfig):
    yield common_nconfig

@pytest.fixture
def common_ndecn_mate_integer(common_nconfig):
    yield common_nconfig

@pytest.fixture
def common_ndecn_mate_real(common_nconfig):
    yield common_nconfig

@pytest.fixture
def common_decn_space_lower_mate_subset(common_nconfig, common_ndecn_mate_subset):
    out = numpy.repeat(0, common_ndecn_mate_subset)
    yield out

@pytest.fixture
def common_decn_space_lower_mate_binary(common_nconfig, common_ndecn_mate_binary):
    out = numpy.repeat(0, common_ndecn_mate_binary)
    yield out

@pytest.fixture
def common_decn_space_lower_mate_integer(common_nconfig, common_ndecn_mate_integer):
    out = numpy.repeat(0, common_ndecn_mate_integer)
    yield out

@pytest.fixture
def common_decn_space_lower_mate_real(common_nconfig, common_ndecn_mate_real):
    out = numpy.repeat(0, common_ndecn_mate_real)
    yield out

@pytest.fixture
def common_decn_space_upper_mate_subset(common_nconfig, common_ndecn_mate_subset):
    out = numpy.repeat(common_nconfig-1, common_ndecn_mate_subset)
    yield out

@pytest.fixture
def common_decn_space_upper_mate_binary(common_nconfig, common_ndecn_mate_binary):
    out = numpy.repeat(1, common_ndecn_mate_binary)
    yield out

@pytest.fixture
def common_decn_space_upper_mate_integer(common_nconfig, common_ndecn_mate_integer):
    out = numpy.repeat(common_nconfig-1, common_ndecn_mate_integer)
    yield out

@pytest.fixture
def common_decn_space_upper_mate_real(common_nconfig, common_ndecn_mate_real):
    out = numpy.repeat(common_nconfig-1, common_ndecn_mate_real)
    yield out

@pytest.fixture
def common_decn_space_mate_subset(common_nconfig):
    out = numpy.arange(common_nconfig)
    yield out

@pytest.fixture
def common_decn_space_mate_binary(common_decn_space_lower_mate_binary, common_decn_space_upper_mate_binary):
    out = numpy.stack([common_decn_space_lower_mate_binary,common_decn_space_upper_mate_binary])
    yield out

@pytest.fixture
def common_decn_space_mate_integer(common_decn_space_lower_mate_integer, common_decn_space_upper_mate_integer):
    out = numpy.stack([common_decn_space_lower_mate_integer,common_decn_space_upper_mate_integer])
    yield out

@pytest.fixture
def common_decn_space_mate_real(common_decn_space_lower_mate_real, common_decn_space_upper_mate_real):
    out = numpy.stack([common_decn_space_lower_mate_real,common_decn_space_upper_mate_real])
    yield out

@pytest.fixture
def common_decn_space_xmap(common_nconfig, common_ntaxa, common_nparent):
    out = numpy.random.randint(0, common_ntaxa, (common_nconfig, common_nparent))
    yield out

@pytest.fixture
def common_nobj():
    yield 1

@pytest.fixture
def common_obj_wt():
    yield None

@pytest.fixture
def common_nineqcv():
    yield None

@pytest.fixture
def common_ineqcv_wt():
    yield None

@pytest.fixture
def common_neqcv():
    yield None

@pytest.fixture
def common_eqcv_wt():
    yield None

@pytest.fixture
def common_nsoln():
    yield 1

@pytest.fixture
def common_soln_decn_mate_subset(common_nconfig, common_ndecn_mate_subset, common_nsoln):
    out = numpy.random.randint(0, common_nconfig, (common_nsoln, common_ndecn_mate_subset))
    yield out

@pytest.fixture
def common_soln_decn_mate_binary(common_nconfig, common_ndecn_mate_binary, common_nsoln):
    out = numpy.random.randint(0, 2, (common_nsoln, common_ndecn_mate_binary))
    yield out

@pytest.fixture
def common_soln_decn_mate_integer(common_nconfig, common_ndecn_mate_integer, common_nsoln):
    out = numpy.random.randint(0, common_nconfig, (common_nsoln, common_ndecn_mate_integer))
    yield out

@pytest.fixture
def common_soln_decn_mate_real(common_nconfig, common_ndecn_mate_real, common_nsoln):
    out = numpy.random.random((common_nsoln, common_ndecn_mate_real))
    yield out

@pytest.fixture
def common_soln_obj(common_nsoln, common_nobj):
    out = numpy.random.random((common_nsoln, common_nobj))
    yield out

@pytest.fixture
def common_soln_ineqcv():
    yield None

@pytest.fixture
def common_soln_eqcv():
    yield None


