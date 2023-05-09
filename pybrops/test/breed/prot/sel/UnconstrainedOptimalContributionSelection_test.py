import os
import numpy
import pytest
from numpy.random import Generator
from numpy.random import PCG64
from matplotlib import pyplot
from matplotlib import animation

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_abstract_function
from pybrops.test import assert_abstract_property
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.opt.algo.NSGA3UnityConstraintGeneticAlgorithm import NSGA3UnityConstraintGeneticAlgorithm
from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.breed.prot.sel.transfn import trans_sum
from pybrops.breed.prot.sel.transfn import trans_ndpt_to_vec_dist
from pybrops.breed.prot.sel.UnconstrainedOptimalContributionSelection import OptimalContributionSelection
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.bvmat.DenseEstimatedBreedingValueMatrix import DenseEstimatedBreedingValueMatrix
from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import DenseMolecularCoancestryMatrix
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
######################## Genotypes #########################
############################################################
@pytest.fixture
def mat_int8():
    yield numpy.int8([
       [[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]],

       [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0]]
    ])

@pytest.fixture
def mat_chrgrp():
    yield numpy.int64([
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2
    ])

@pytest.fixture
def mat_phypos():
    yield numpy.int64([
         1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
        11, 12, 13, 14, 15, 16, 17, 18, 19, 20
    ])

@pytest.fixture
def mat_genpos():
    yield numpy.float64([
        0.13, 0.32, 0.53, 0.54, 0.55, 0.61, 0.63, 0.7 , 0.75, 0.96,
        0.14, 0.16, 0.26, 0.31, 0.31, 0.68, 0.7 , 0.74, 0.75, 0.91
    ])

@pytest.fixture
def mat_taxa():
    yield numpy.object_([
        'Line01', 'Line02', 'Line03', 'Line04', 'Line05',
        'Line06', 'Line07', 'Line08', 'Line09', 'Line10',
        'Line11', 'Line12', 'Line13', 'Line14', 'Line15',
        'Line16', 'Line17', 'Line18', 'Line19', 'Line20'
    ])

@pytest.fixture
def mat_taxa_grp():
    yield numpy.int64([
        1, 1, 1, 1, 1,
        2, 2, 2, 2, 2,
        3, 3, 3, 3, 3,
        4, 4, 4, 4, 4
    ])

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
###################### Genomic model #######################
############################################################
@pytest.fixture
def mat_beta():
    yield numpy.float64([
        [25.6, 13.4]
    ])
    # yield numpy.float64([[1.4, 2.5, 7.2]])

@pytest.fixture
def mat_u_misc():
    yield None

@pytest.fixture
def mat_u_a():
    yield numpy.float64([
        [ 1.25, -0.68],
        [-0.02, -1.09],
        [ 0.21, -0.5 ],
        [-2.84,  0.64],
        [-1.37, -0.81],
        [-2.06,  2.22],
        [ 1.52, -0.21],
        [-0.23, -1.78],
        [ 1.04, -0.55],
        [-0.77, -1.4 ],
        [-0.44,  0.89],
        [ 0.12, -0.87],
        [-0.44, -0.55],
        [ 1.36,  0.73],
        [ 1.04,  1.22],
        [-0.05,  0.82],
        [ 0.93,  0.73],
        [-0.89,  1.21],
        [ 0.05, -1.19],
        [-1.27, -2.  ]
    ])

@pytest.fixture
def trait():
    yield numpy.object_(["protein", "yield"])

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
################## Breeding values model ###################
############################################################
@pytest.fixture
def bvmat(dalgmod, dgmat):
    yield dalgmod.gebv(dgmat)

############################################################
############### OptimalContributionSelection ###############
############################################################
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
def inbfn():
    def fn(t_cur, t_max):
        return 0.75
    yield fn

@pytest.fixture
def cmatcls():
    yield DenseMolecularCoancestryMatrix

@pytest.fixture
def method():
    yield "single"

@pytest.fixture
def ocs(nparent, ncross, nprogeny, inbfn, cmatcls, method, rng):
    yield OptimalContributionSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        inbfn = inbfn,
        cmatcls = cmatcls,
        method = method,
        objfn_trans = trans_sum, # sum of two traits
        objfn_trans_kwargs = None, # no kwargs
        objfn_wt = 1.0, # maximizing
        rng = rng
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(OptimalContributionSelection)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(OptimalContributionSelection, "__init__")

def test_select_is_concrete():
    assert_concrete_method(OptimalContributionSelection, "select")

def test_objfn_is_concrete():
    assert_concrete_method(OptimalContributionSelection, "objfn")

def test_objfn_vec_is_concrete():
    assert_concrete_method(OptimalContributionSelection, "objfn_vec")

def test_pareto_is_concrete():
    assert_concrete_method(OptimalContributionSelection, "pareto")

def test_objfn_static_is_concrete():
    assert_concrete_method(OptimalContributionSelection, "objfn_static")

def test_objfn_vec_static_is_concrete():
    assert_concrete_method(OptimalContributionSelection, "objfn_vec_static")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_select_single(nparent, ncross, nprogeny, inbfn, cmatcls, method, rng, dpgmat, dgmat, bvmat, dalgmod):
    # make selection object
    sel = OptimalContributionSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        inbfn = inbfn,
        cmatcls = cmatcls,
        method = "single",
        objfn_trans = trans_sum, # sum of two traits
        objfn_trans_kwargs = {}, # no kwargs
        objfn_wt = 1.0, # maximizing
        rng = rng
    )
    # make selections
    pgmat, sel, ncross, nprogeny = sel.select(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )

    assert sel.ndim == 1

def test_select_single_objfn_trans_RuntimeError(nparent, ncross, nprogeny, inbfn, cmatcls, method, rng, dpgmat, dgmat, bvmat, dalgmod):
    # make selection object
    sel = OptimalContributionSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        inbfn = inbfn,
        cmatcls = cmatcls,
        method = "single",
        objfn_trans = None, # sum of two traits
        objfn_trans_kwargs = None, # no kwargs
        objfn_wt = 1.0, # maximizing function
        ndset_trans = trans_ndpt_to_vec_dist,
        ndset_trans_kwargs = {
            "objfn_wt": numpy.array([1.0, 1.0, 1.0]),   # all objectives maximizing
            "wt": numpy.array([1./3., 1./3., 1./3.])    # 1/3 equal weight to all
        },
        rng = rng
    )

    # make sure this raises an error due to lack of shape compatibility
    with pytest.raises(RuntimeError):
        # make selections
        pgmat, sel, ncross, nprogeny = sel.select(
            pgmat = dpgmat,
            gmat = dgmat,
            ptdf = None,
            bvmat = bvmat,
            gpmod = dalgmod,
            t_cur = 0,
            t_max = 20
        )

def test_select_pareto(nparent, ncross, nprogeny, inbfn, cmatcls, method, rng, dpgmat, dgmat, bvmat, dalgmod):
    # make selection object
    sel = OptimalContributionSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        inbfn = inbfn,
        cmatcls = cmatcls,
        method = "pareto",
        objfn_trans = None, # sum of two traits
        objfn_trans_kwargs = None, # no kwargs
        objfn_wt = numpy.array([1.0, 1.0, 1.0]), # maximizing function
        ndset_trans = trans_ndpt_to_vec_dist,
        ndset_trans_kwargs = {
            "objfn_wt": numpy.array([1.0, 1.0, 1.0]),   # all objectives maximizing
            "wt": numpy.array([1./3., 1./3., 1./3.])    # 1/3 equal weight to all
        },
        rng = rng
    )
    # make selections
    pgmat, sel, ncross, nprogeny = sel.select(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )

    assert sel.ndim == 1

def test_objfn_is_function(ocs, dgmat, bvmat, dalgmod):
    objfn = ocs.objfn(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )

    assert callable(objfn)

def test_objfn_multiobjective(nparent, ncross, nprogeny, inbfn, cmatcls, method, rng, dpgmat, dgmat, bvmat, dalgmod):
    # make selection object
    sel = OptimalContributionSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        inbfn = inbfn,
        cmatcls = cmatcls,
        method = "pareto",
        objfn_trans = None, # sum of two traits
        objfn_trans_kwargs = None, # no kwargs
        objfn_wt = numpy.array([1.0, 1.0, 1.0]), # maximizing function
        ndset_trans = trans_ndpt_to_vec_dist,
        ndset_trans_kwargs = {
            "objfn_wt": numpy.array([1.0, 1.0, 1.0]),   # all objectives maximizing
            "wt": numpy.array([1./3., 1./3., 1./3.])    # 1/3 equal weight to all
        },
        rng = rng
    )

    # make objective function
    objfn = sel.objfn(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )

    # get number of taxa
    ntaxa = dgmat.ntaxa

    # make contribution selection
    x = numpy.random.uniform(0,1,ntaxa)
    x *= (1.0/x.sum())

    out = objfn(x)

    assert out.ndim == 1
    assert out.shape[0] > 1

def test_objfn_vec_is_callable(ocs, dgmat, bvmat, dalgmod):
    objfn_vec = ocs.objfn_vec(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )

    assert callable(objfn_vec)

def test_objfn_vec_multiobjective(nparent, ncross, nprogeny, inbfn, cmatcls, method, rng, dpgmat, dgmat, bvmat, dalgmod):
    # make selection object
    sel = OptimalContributionSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        inbfn = inbfn,
        cmatcls = cmatcls,
        method = "pareto",
        objfn_trans = None, # sum of two traits
        objfn_trans_kwargs = None, # no kwargs
        objfn_wt = numpy.array([1.0, 1.0, 1.0]), # maximizing function
        ndset_trans = trans_ndpt_to_vec_dist,
        ndset_trans_kwargs = {
            "objfn_wt": numpy.array([1.0, 1.0, 1.0]),   # all objectives maximizing
            "wt": numpy.array([1./3., 1./3., 1./3.])    # 1/3 equal weight to all
        },
        rng = rng
    )

    # make objective function
    objfn = sel.objfn_vec(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )

    # get number of taxa
    ntaxa = dgmat.ntaxa

    # number of simulated solutions
    nsoln = 5

    # make contribution selection
    x = numpy.random.uniform(0,1,(nsoln,ntaxa))
    x *= (1.0/x.sum(1))[:,None]

    out = objfn(x)

    assert out.ndim == 2
    assert out.shape[0] == nsoln
    assert out.shape[1] > 1

def test_pareto(nparent, ncross, nprogeny, inbfn, cmatcls, method, rng, dpgmat, dgmat, bvmat, dalgmod):
    # make selection object
    sel = OptimalContributionSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        inbfn = inbfn,
        cmatcls = cmatcls,
        method = "pareto",
        objfn_trans = None, # sum of two traits
        objfn_trans_kwargs = None, # no kwargs
        objfn_wt = numpy.array([1.0, 1.0, 1.0]), # maximizing function
        ndset_trans = trans_ndpt_to_vec_dist,
        ndset_trans_kwargs = {
            "objfn_wt": numpy.array([1.0, 1.0, 1.0]),   # all objectives maximizing
            "wt": numpy.array([1./3., 1./3., 1./3.])    # 1/3 equal weight to all
        },
        moalgo = NSGA3UnityConstraintGeneticAlgorithm(
            ngen = 1000,            # number of generations to evolve
            mu = 100,               # number of parents in population
            lamb = 100,             # number of progeny to produce
            cxeta = 30.0,           # crossover variance parameter
            muteta = 20.0,          # mutation crossover parameter
            refpnts = None,         # hyperplane reference points
            save_logbook = False,   # whether to save logs or not
            rng = rng               # PRNG source
        ),
        rng = rng
    )

    # get pareto frontier
    frontier, sel_config = sel.pareto(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )

    xdata = frontier[:,0]
    ydata = frontier[:,1]
    zdata = frontier[:,2]

    xlabel = "inbreeding"
    ylabel = dalgmod.trait[0]
    zlabel = dalgmod.trait[1]

    # create static figure
    fig = pyplot.figure()
    ax = pyplot.axes(projection = '3d')
    ax.scatter3D(xdata, ydata, zdata)
    ax.set_title("Optimal Contribution Selection Test Pareto Frontier")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    pyplot.savefig("OCS_3d_frontier.png", dpi = 250)

    # create animation
    fig = pyplot.figure()
    ax = pyplot.axes(projection = '3d')

    def init():
        ax.scatter3D(xdata, ydata, zdata)
        ax.set_title("Optimal Contribution Selection Test Pareto Frontier")
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
        s = outdir + "/" + "OCS_3d_frontier_" + str(i).zfill(3) + ".png"
        pyplot.savefig(s, dpi = 250)

    # does not want to work!
    # anim = animation.FuncAnimation(fig, animate, init_func = init, frames = 100, interval = 100, blit = True)
    # writer = animation.PillowWriter(fps=30)
    # writer = animation.FFMpegWriter(fps = 10)
    # anim.save('OCS_3d_frontier.gif', writer = writer)


# def test_objfn_static_multiobjective(ocs, mat_int8, mat_u_a):
#     Z = mat_int8.sum(0) # (n,p) genotypes {0,1,2}
#     u = mat_u_a         # (p,t) regression coefficients
#     Y = Z@u             # (n,t) values
#
#     for i,taxon_bv in enumerate(Y):
#         numpy.testing.assert_almost_equal(
#             taxon_bv,
#             ocs.objfn_static([i], Z, u, None, None)
#         )
#
# def test_objfn_vec_static_multiobjective(ocs, dgmat, bvmat, dalgmod, mat_int8, mat_u_a):
#     Z = mat_int8.sum(0) # (n,p) genotypes {0,1,2}
#     u = mat_u_a         # (p,t) regression coefficients
#     Y = Z@u             # (n,t) values
#
#     for i,taxon_bv in enumerate(Y):
#         numpy.testing.assert_almost_equal(
#             numpy.stack([taxon_bv, taxon_bv]),
#             ocs.objfn_vec_static([[i],[i]], Z, u, None, None)
#         )
