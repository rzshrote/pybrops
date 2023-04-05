from matplotlib import pyplot
import numpy
import pytest
from numpy.random import Generator
from numpy.random import PCG64
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.model.vmat.fcty.DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory import DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_abstract_function
from pybrops.test import assert_abstract_property
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.breed.prot.sel.transfn import trans_ndpt_to_vec_dist, trans_sum
from pybrops.breed.prot.sel.UsefulnessCriterionSelection import UsefulnessCriterionSelection
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
######################## Genotypes #########################
############################################################
@pytest.fixture
def mat_ntaxa():
    # must be divisible by 4
    yield 40

@pytest.fixture
def mat_nvrnt():
    # must be divisible by 2
    yield 100

@pytest.fixture
def mat_int8(mat_ntaxa, mat_nvrnt):
    yield numpy.random.binomial(1,0.1,(2,mat_ntaxa,mat_nvrnt)).astype('int8')

@pytest.fixture
def mat_chrgrp(mat_nvrnt):
    yield numpy.repeat([1,2], mat_nvrnt // 2)

@pytest.fixture
def mat_phypos(mat_nvrnt):
    yield numpy.arange(1, mat_nvrnt+1)

@pytest.fixture
def mat_genpos(mat_nvrnt):
    out = numpy.random.random(mat_nvrnt)
    out.sort()
    yield out

@pytest.fixture
def mat_taxa(mat_ntaxa):
    yield numpy.array(["Line" + str(i).zfill(3) for i in range(1,mat_ntaxa+1)], dtype = 'object')

@pytest.fixture
def mat_taxa_grp(mat_ntaxa):
    yield numpy.repeat([1,2,3,4], mat_ntaxa // 4)

@pytest.fixture
def dpgmat(mat_int8, mat_chrgrp, mat_phypos, mat_genpos, mat_taxa, mat_taxa_grp):
    out = DensePhasedGenotypeMatrix(
        mat = mat_int8,
        vrnt_chrgrp = mat_chrgrp,
        vrnt_phypos = mat_phypos,
        vrnt_genpos = mat_genpos,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp
    )
    out.group()
    yield out

@pytest.fixture
def dgmat(dpgmat):
    dugt = DenseUnphasedGenotyping()
    out = dugt.genotype(dpgmat)
    yield out

############################################################
###################### Genomic Models ######################
############################################################
@pytest.fixture
def mat_ntrait():
    yield 2

@pytest.fixture
def mat_beta(mat_ntrait):
    yield numpy.random.uniform(10, 50, (1,mat_ntrait))

@pytest.fixture
def mat_u_misc():
    yield None

@pytest.fixture
def mat_u_a(mat_nvrnt, mat_ntrait):
    yield numpy.random.normal(0,1,(mat_nvrnt,mat_ntrait))

@pytest.fixture
def trait(mat_ntrait):
    yield numpy.array(["Trait" + str(i).zfill(2) for i in range(1,mat_ntrait+1)], dtype = 'object')

@pytest.fixture
def model_name():
    yield "test_dalgmod"

@pytest.fixture
def params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def dalgmod(mat_beta, mat_u_misc, mat_u_a, trait, model_name, params):
    yield DenseAdditiveLinearGenomicModel(
        beta = mat_beta,
        u_misc = mat_u_misc,
        u_a = mat_u_a,
        trait = trait,
        model_name = model_name,
        params = params
    )

############################################################
############### UsefulnessCriterionSelection ###############
############################################################
@pytest.fixture
def nconfig():
    yield 5

@pytest.fixture
def nparent():
    yield 2

@pytest.fixture
def ncross():
    yield 1

@pytest.fixture
def nprogeny():
    yield 10

@pytest.fixture
def rng():
    yield Generator(PCG64(192837465))

@pytest.fixture
def method():
    yield "single"

@pytest.fixture
def ucs(nconfig, nparent, ncross, nprogeny, method, rng):
    yield UsefulnessCriterionSelection(
        nconfig = nconfig,
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        nself = 0,
        upper_percentile = 0.1,
        vmatfcty = DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory(),
        gmapfn = HaldaneMapFunction(),
        method = method,
        objfn_trans = trans_sum, # sum of two traits
        objfn_trans_kwargs = {}, # no kwargs
        objfn_wt = 1.0, # maximizing
        rng = rng
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(UsefulnessCriterionSelection)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(UsefulnessCriterionSelection, "__init__")

def test_select_is_concrete():
    assert_concrete_method(UsefulnessCriterionSelection, "select")

def test_objfn_is_concrete():
    assert_concrete_method(UsefulnessCriterionSelection, "objfn")

def test_objfn_vec_is_concrete():
    assert_concrete_method(UsefulnessCriterionSelection, "objfn_vec")

def test_pareto_is_concrete():
    assert_concrete_method(UsefulnessCriterionSelection, "pareto")

def test_objfn_static_is_concrete():
    assert_concrete_method(UsefulnessCriterionSelection, "objfn_static")

def test_objfn_vec_static_is_concrete():
    assert_concrete_method(UsefulnessCriterionSelection, "objfn_vec_static")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################
def test_calc_xmap(ucs, mat_ntaxa, nparent):
    xmap = ucs._calc_xmap(mat_ntaxa)
    if ucs.unique_parents:
        assert xmap.shape[0] == (mat_ntaxa * (mat_ntaxa-1))//2 # assumes 2 parnets
    else:
        assert xmap.shape[0] == (mat_ntaxa * (mat_ntaxa+1))//2
    assert xmap.shape[1] == nparent

def test_calc_uc(ucs, dpgmat, dalgmod, mat_ntaxa, mat_ntrait):
    xmap = ucs._calc_xmap(mat_ntaxa)
    uc = ucs._calc_uc(dpgmat, dalgmod, xmap)
    if ucs.unique_parents:
        assert uc.shape[0] == (mat_ntaxa * (mat_ntaxa-1))//2 # assumes 2 parnets
    else:
        assert uc.shape[0] == (mat_ntaxa * (mat_ntaxa+1))//2
    assert uc.shape[1] == mat_ntrait

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_select_single(ucs, ncross, nprogeny, dpgmat, dgmat, dalgmod):
    # make selections
    sel_pgmat, sel, sel_ncross, sel_nprogeny = ucs.select(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )
    assert id(sel_pgmat) == id(dpgmat)
    assert sel.ndim == 1
    assert sel_ncross == ncross
    assert sel_nprogeny == nprogeny

def test_select_pareto_TypeError(ucs, dpgmat, dalgmod):
    # set to pareto selection
    ucs.method = "pareto"

    # make sure pareto raises errors
    with pytest.raises(TypeError):
        # make selections
        pgmat, sel, ncross, nprogeny = ucs.select(
            pgmat = None,
            gmat = None,
            ptdf = None,
            bvmat = None,
            gpmod = dalgmod,
            t_cur = 0,
            t_max = 20
        )

    # make sure pareto raises errors
    with pytest.raises(TypeError):
        # make selections
        pgmat, sel, ncross, nprogeny = ucs.select(
            pgmat = dpgmat,
            gmat = None,
            ptdf = None,
            bvmat = None,
            gpmod = None,
            t_cur = 0,
            t_max = 20
        )

def test_select_pareto(nconfig, nparent, ncross, nprogeny, rng, dpgmat, dalgmod):
    ucs = UsefulnessCriterionSelection(
        nconfig = nconfig,
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        nself = 0,
        upper_percentile = 0.1,
        vmatfcty = DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory(),
        gmapfn = HaldaneMapFunction(),
        method = "pareto",
        objfn_trans = None, # sum of two traits
        objfn_trans_kwargs = {}, # no kwargs
        objfn_wt = numpy.array([1.0,1.0]), # maximizing
        ndset_trans = trans_ndpt_to_vec_dist,
        ndset_trans_kwargs = {"objfn_wt": numpy.array([1.0,1.0]), "wt": numpy.array([0.5,0.5])},
        ndset_wt = -1.0,
        rng = rng
    )

    # make selections
    sel_pgmat, sel, sel_ncross, sel_nprogeny = ucs.select(
        pgmat = dpgmat,
        gmat = None,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )
    assert id(sel_pgmat) == id(dpgmat)
    assert sel.ndim == 1
    assert sel_ncross == ncross
    assert sel_nprogeny == nprogeny

def test_pareto(nconfig, nparent, ncross, nprogeny, rng, dpgmat, dalgmod):
    ucs = UsefulnessCriterionSelection(
        nconfig = nconfig,
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        nself = 0,
        upper_percentile = 0.1,
        vmatfcty = DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory(),
        gmapfn = HaldaneMapFunction(),
        method = "pareto",
        objfn_trans = None, # sum of two traits
        objfn_trans_kwargs = {}, # no kwargs
        objfn_wt = numpy.array([1.0,1.0]), # maximizing
        rng = rng
    )

    frontier, sel_config = ucs.pareto(
        pgmat = dpgmat,
        gmat = None,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )

    xdata = frontier[:,0]
    ydata = frontier[:,1]
    # zdata = frontier[:,2]

    xlabel = dalgmod.trait[0]
    ylabel = dalgmod.trait[1]

    fig = pyplot.figure()
    ax = pyplot.axes()
    ax.scatter(xdata, ydata)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Usefulness Criterion Selection Test Pareto Frontier")
    # ax = pyplot.axes(projection='3d')
    # ax.scatter3D(xdata, ydata, zdata)
    pyplot.savefig("UCS_2d_frontier.png", dpi = 250)

def test_objfn_is_function(ucs, dpgmat, dalgmod):
    # make objective function
    objfn = ucs.objfn(
        pgmat = dpgmat,
        gmat = None,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )
    assert callable(objfn)

def test_objfn_is_function(ucs, dpgmat, dalgmod):
    # make objective function
    objfn = ucs.objfn_vec(
        pgmat = dpgmat,
        gmat = None,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )
    assert callable(objfn)

def test_objfn_is_multiobjective(ucs, dpgmat, dalgmod, mat_ntaxa, mat_ntrait):
    # set transformation function to None
    ucs.objfn_trans = None
    # make objective function
    objfn = ucs.objfn(
        pgmat = dpgmat,
        gmat = None,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )
    # make random selection
    x = numpy.random.choice(mat_ntaxa, 5)
    out = objfn(x)
    assert out.ndim == 1
    assert out.shape[0] == mat_ntrait

def test_objfn_vec_is_multiobjective(ucs, dpgmat, dalgmod, mat_ntaxa, mat_ntrait):
    # set transformation function to None
    ucs.objfn_trans = None
    # make objective function
    objfn = ucs.objfn_vec(
        pgmat = dpgmat,
        gmat = None,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )
    # make random selection: 7 sets of 5 config
    x = numpy.random.choice(mat_ntaxa, (7,5))
    out = objfn(x)
    assert out.ndim == 2
    assert out.shape[0] == 7
    assert out.shape[1] == mat_ntrait
