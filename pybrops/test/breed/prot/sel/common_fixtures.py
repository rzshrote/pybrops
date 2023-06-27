import random
import numpy
import pytest
from numpy.random import Generator, PCG64
from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.model.vmat.fcty.DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory import DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.breed.prot.sel.prob.trans import trans_ndpt_to_vec_dist


################################ Shape fixtures ################################

@pytest.fixture
def common_ntaxa():
    yield 100 # must be divisible by 4

@pytest.fixture
def common_nvrnt():
    yield 1000

@pytest.fixture
def common_ntrait():
    yield 2

@pytest.fixture
def common_nself():
    yield 0

@pytest.fixture
def common_unique_parents():
    yield True

@pytest.fixture
def common_nchr():
    yield 2

@pytest.fixture
def common_ngrp():
    yield 4

@pytest.fixture
def common_upper_percentile():
    yield 0.1

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

########################### Genotype matrix fixtures ###########################

@pytest.fixture
def common_gmat(common_pgmat):
    dugt = DenseUnphasedGenotyping()
    out = dugt.genotype(common_pgmat)
    yield out

@pytest.fixture
def common_algmod_u_a_mean(common_ntrait):
    out = numpy.repeat(0.0, common_ntrait)
    yield out

@pytest.fixture
def common_algmod_u_a_cov(common_ntrait):
    pass

@pytest.fixture
def common_algmod_beta(common_ntrait):
    out = numpy.random.uniform(10, 30, (1, common_ntrait))
    yield out

########################## Variance matrix fixtures ###########################
@pytest.fixture
def common_vmatfcty():
    out = DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory()
    yield out

############################# Genetic map fixtures #############################
@pytest.fixture
def common_gmapfn():
    out = HaldaneMapFunction()
    yield out

############################ Genomic model fixtures ############################
@pytest.fixture
def common_algmod(common_nvrnt, common_ntrait):
    beta = numpy.random.uniform(10, 30, (1,common_ntrait))
    u_misc = None
    u_a = numpy.random.normal(0, 1, (common_nvrnt,common_ntrait))
    trait = numpy.array(["Trait"+str(i).zfill(2) for i in range(1,common_ntrait+1)], dtype=object)
    model_name = "test_dalgmod"
    params = {"a" : 0, "b" : 1}
    out = DenseAdditiveLinearGenomicModel(
        beta = beta,
        u_misc = u_misc,
        u_a = u_a,
        trait = trait,
        model_name = model_name,
        params = params
    )
    yield out

@pytest.fixture
def common_gpmod(common_algmod):
    yield common_algmod

######################### Phenotype dataframe fixtures #########################
@pytest.fixture
def common_ptdf():
    yield None

########################### Breeding value fixtures ############################
@pytest.fixture
def common_bvmat(common_algmod, common_gmat):
    yield common_algmod.gebv(common_gmat)

@pytest.fixture
def common_t_cur():
    yield 0

@pytest.fixture
def common_t_max():
    yield 20

########################## SelectionProblem fixtures ###########################
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

@pytest.fixture
def common_nobj():
    yield 2

@pytest.fixture
def common_obj_wt():
    yield numpy.array([1,-1], dtype=float)

@pytest.fixture
def common_obj_trans():
    yield None

@pytest.fixture
def common_obj_trans_kwargs():
    yield None

@pytest.fixture
def common_nineqcv():
    yield 0

@pytest.fixture
def common_ineqcv_wt():
    yield numpy.array([], dtype=float)

@pytest.fixture
def common_ineqcv_trans():
    yield None

@pytest.fixture
def common_ineqcv_trans_kwargs():
    yield None

@pytest.fixture
def common_neqcv():
    yield 0

@pytest.fixture
def common_eqcv_wt():
    yield numpy.array([], dtype=float)

@pytest.fixture
def common_eqcv_trans():
    yield None

@pytest.fixture
def common_eqcv_trans_kwargs():
    yield None

@pytest.fixture
def common_ndset_wt():
    yield -1.0

@pytest.fixture
def common_ndset_trans():
    yield trans_ndpt_to_vec_dist

@pytest.fixture
def common_ndset_trans_kwargs(common_obj_wt):
    yield {"obj_wt": common_obj_wt, "vec_wt": numpy.array([0.5,0.5])}

@pytest.fixture
def common_rng():
    yield Generator(PCG64(192837465))

@pytest.fixture
def common_soalgo():
    yield None

@pytest.fixture
def common_moalgo():
    yield None
