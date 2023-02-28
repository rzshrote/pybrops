import os
from matplotlib import pyplot
import numpy
import pytest
from numpy.random import Generator
from numpy.random import PCG64
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.model.vmat.fcty.DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory import DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory
from pybrops.popgen.cmat.fcty.DenseMolecularCoancestryMatrixFactory import DenseMolecularCoancestryMatrixFactory
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction

from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.breed.prot.sel.transfn import trans_ndpt_to_vec_dist, trans_sum, trans_sum_inbmax_penalty
from pybrops.breed.prot.sel.BinaryOptimalContributionSelection import BinaryOptimalContributionSelection
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

@pytest.fixture
def bvmat(dalgmod, dpgmat):
    yield dalgmod.gebv(dpgmat)

############################################################
############### BinaryOptimalContributionSelection ###############
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
def inbfn():
    def fn(t_cur, t_max):
        return 1 - (t_cur/t_max)
    yield fn

@pytest.fixture
def cmatfcty():
    yield DenseMolecularCoancestryMatrixFactory()

@pytest.fixture
def bocs(nparent, ncross, nprogeny, inbfn, cmatfcty, method, rng):
    yield BinaryOptimalContributionSelection(
        nparent = nparent, 
        ncross = ncross, 
        nprogeny = nprogeny,
        inbfn = inbfn,
        cmatfcty = cmatfcty,
        method = method,
        objfn_trans = trans_sum_inbmax_penalty, 
        objfn_trans_kwargs = {"penalty_wt": -10}, 
        objfn_wt = 1.0,
        ndset_trans = None, 
        ndset_trans_kwargs = None, 
        ndset_wt = -1.0,
        rng = rng, 
        soalgo = None,
        moalgo = None
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(BinaryOptimalContributionSelection)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(BinaryOptimalContributionSelection, "__init__")

def test_select_is_concrete():
    generic_assert_concrete_method(BinaryOptimalContributionSelection, "select")

def test_objfn_is_concrete():
    generic_assert_concrete_method(BinaryOptimalContributionSelection, "objfn")

def test_objfn_vec_is_concrete():
    generic_assert_concrete_method(BinaryOptimalContributionSelection, "objfn_vec")

def test_pareto_is_concrete():
    generic_assert_concrete_method(BinaryOptimalContributionSelection, "pareto")

def test_objfn_static_is_concrete():
    generic_assert_concrete_method(BinaryOptimalContributionSelection, "objfn_static")

def test_objfn_vec_static_is_concrete():
    generic_assert_concrete_method(BinaryOptimalContributionSelection, "objfn_vec_static")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_select_single(bocs, ncross, nprogeny, dpgmat, dgmat, dalgmod, bvmat):
    # make selections
    sel_pgmat, sel, sel_ncross, sel_nprogeny = bocs.select(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    assert id(sel_pgmat) == id(dpgmat)
    assert sel.ndim == 1
    assert sel_ncross == ncross
    assert sel_nprogeny == nprogeny

def test_select_pareto_TypeError(bocs, dpgmat, dgmat, bvmat):
    # set to pareto selection
    bocs.method = "pareto"

    # make sure pareto raises errors
    with pytest.raises(TypeError):
        # make selections
        pgmat, sel, ncross, nprogeny = bocs.select(
            pgmat = None,
            gmat = dgmat,
            ptdf = None,
            bvmat = bvmat,
            gpmod = None,
            t_cur = 0,
            t_max = 20
        )

    # make sure pareto raises errors
    with pytest.raises(TypeError):
        # make selections
        pgmat, sel, ncross, nprogeny = bocs.select(
            pgmat = dpgmat,
            gmat = None,
            ptdf = None,
            bvmat = bvmat,
            gpmod = None,
            t_cur = 0,
            t_max = 20
        )

    # make sure pareto raises errors
    with pytest.raises(TypeError):
        # make selections
        pgmat, sel, ncross, nprogeny = bocs.select(
            pgmat = dpgmat,
            gmat = dgmat,
            ptdf = None,
            bvmat = None,
            gpmod = None,
            t_cur = 0,
            t_max = 20
        )

def test_select_pareto(nparent, ncross, nprogeny, inbfn, cmatfcty, rng, dpgmat, dgmat, bvmat):
    bocs = BinaryOptimalContributionSelection(
        nparent = nparent, 
        ncross = ncross, 
        nprogeny = nprogeny,
        inbfn = inbfn,
        cmatfcty = cmatfcty,
        method = "pareto",
        objfn_trans = None, 
        objfn_trans_kwargs = None, 
        objfn_wt = numpy.array([-1.0,1.0,1.0]), # min inbreeding, max bv
        ndset_trans = trans_ndpt_to_vec_dist,
        ndset_trans_kwargs = {"objfn_wt": numpy.array([-1.0,1.0,1.0]), "wt": numpy.array([0.33,0.33,0.33])},
        ndset_wt = -1.0,
        rng = rng
    )

    # make selections
    sel_pgmat, sel, sel_ncross, sel_nprogeny = bocs.select(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    assert id(sel_pgmat) == id(dpgmat)
    assert sel.ndim == 1
    assert sel_ncross == ncross
    assert sel_nprogeny == nprogeny

def test_pareto(nconfig, nparent, ncross, nprogeny, rng, dpgmat, dgmat, bvmat, inbfn, cmatfcty):
    bocs = BinaryOptimalContributionSelection(
        nparent = nparent, 
        ncross = ncross, 
        nprogeny = nprogeny,
        inbfn = inbfn,
        cmatfcty = cmatfcty,
        method = "pareto",
        objfn_trans = None, 
        objfn_trans_kwargs = None, 
        objfn_wt = numpy.array([-1.0,1.0,1.0]), # min inbreeding, max bv
        ndset_trans = trans_ndpt_to_vec_dist,
        ndset_trans_kwargs = {"objfn_wt": numpy.array([-1.0,1.0,1.0]), "wt": numpy.array([0.33,0.33,0.33])},
        ndset_wt = -1.0,
        rng = rng
    )

    frontier, sel_config = bocs.pareto(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )

    xdata = frontier[:,0]
    ydata = frontier[:,1]
    zdata = frontier[:,2]

    xlabel = "inbreeding"
    ylabel = bvmat.trait[0]
    zlabel = bvmat.trait[1]

    # create static figure
    fig = pyplot.figure()
    ax = pyplot.axes(projection = '3d')
    ax.scatter3D(xdata, ydata, zdata)
    ax.set_title("Binary Optimal Contribution Selection Test Pareto Frontier")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    pyplot.savefig("BOCS_3d_frontier.png", dpi = 250)

    # create animation
    fig = pyplot.figure()
    ax = pyplot.axes(projection = '3d')

    def init():
        ax.scatter3D(xdata, ydata, zdata)
        ax.set_title("Binary Optimal Contribution Selection Test Pareto Frontier")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        return fig,

    def animate(i):
        ax.view_init(elev = 30., azim = 3.6 * i)
        # ax.view_init(elev = 30., azim = i)
        return fig,

    # create and same animation
    outdir = "frames"
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    init()
    for i in range(100):
    # for i in range(360):
        animate(i)
        s = outdir + "/" + "BOCS_3d_frontier_" + str(i).zfill(3) + ".png"
        pyplot.savefig(s, dpi = 250)

def test_objfn_is_function(bocs, dpgmat, dgmat, bvmat):
    # make objective function
    objfn = bocs.objfn(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    assert callable(objfn)

def test_objfn_is_function(bocs, dpgmat, dgmat, bvmat):
    # make objective function
    objfn = bocs.objfn_vec(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    assert callable(objfn)

def test_objfn_is_multiobjective(bocs, dpgmat, dgmat, bvmat, mat_ntaxa, mat_ntrait):
    # set transformation function to None
    bocs.objfn_trans = None
    # make objective function
    objfn = bocs.objfn(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    # make random selection
    x = numpy.random.choice(mat_ntaxa, 5)
    out = objfn(x)
    assert out.ndim == 1
    assert out.shape[0] == (1 + mat_ntrait)

def test_objfn_vec_is_multiobjective(bocs, dpgmat, dgmat, bvmat, mat_ntaxa, mat_ntrait):
    # set transformation function to None
    bocs.objfn_trans = None
    # make objective function
    objfn = bocs.objfn_vec(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    # make random selection: 7 sets of 5 config
    x = numpy.random.choice(mat_ntaxa, (7,5))
    out = objfn(x)
    assert out.ndim == 2
    assert out.shape[0] == 7
    assert out.shape[1] == (1 + mat_ntrait)
