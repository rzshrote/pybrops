import numpy
import pytest

from numpy.random import Generator
from numpy.random import PCG64
from matplotlib import pyplot
from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_abstract_function
from pybrops.test import assert_abstract_property
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.breed.prot.sel.OptimalPopulationValueSelection import OptimalPopulationValueSelection
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.bvmat.DenseEstimatedBreedingValueMatrix import DenseEstimatedBreedingValueMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.opt.algo.SteepestAscentSetHillClimber import SteepestAscentSetHillClimber
from pybrops.breed.prot.sel.transfn import trans_ndpt_to_vec_dist, trans_sum

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
###################### Genomic model #######################
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
################## Breeding values model ###################
############################################################
@pytest.fixture
def bvmat(dalgmod, dpgmat):
    yield dalgmod.gebv(dpgmat)

############################################################
############# OptimalPopulationValueSelection ##############
############################################################
@pytest.fixture
def nparent():
    yield 4

@pytest.fixture
def ncross():
    yield 1

@pytest.fixture
def nprogeny():
    yield 10

@pytest.fixture
def nhaploblk():
    yield 4

@pytest.fixture
def objfn_trans():
    yield trans_sum

@pytest.fixture
def method():
    yield "single"

@pytest.fixture
def objfn_wt():
    yield 1.0

@pytest.fixture
def rng():
    yield Generator(PCG64(549824248))

@pytest.fixture
def opvs(nparent, ncross, nprogeny, nhaploblk, objfn_trans, method, rng):
    yield OptimalPopulationValueSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        nhaploblk = nhaploblk,
        objfn_trans = objfn_trans,
        method = method,
        rng = rng
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(OptimalPopulationValueSelection)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(OptimalPopulationValueSelection, "__init__")

def test_select_is_concrete():
    assert_concrete_method(OptimalPopulationValueSelection, "select")

def test_objfn_is_concrete():
    assert_concrete_method(OptimalPopulationValueSelection, "objfn")

def test_objfn_vec_is_concrete():
    assert_concrete_method(OptimalPopulationValueSelection, "objfn_vec")

def test_pareto_is_concrete():
    assert_concrete_method(OptimalPopulationValueSelection, "pareto")

def test_objfn_static_is_concrete():
    assert_concrete_method(OptimalPopulationValueSelection, "objfn_static")

def test_objfn_vec_static_is_concrete():
    assert_concrete_method(OptimalPopulationValueSelection, "objfn_vec_static")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################
@pytest.fixture
def test_calc_hmat(opvs, dpgmat, dalgmod, mat_ntaxa, nhaploblk, mat_ntrait):
    hmat = opvs._calc_hmat(dpgmat, dalgmod)
    assert hmat.shape[0] == dpgmat.nphase
    assert hmat.shape[1] == mat_ntaxa
    assert hmat.shape[2] == nhaploblk
    assert hmat.shape[3] == mat_ntrait

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_select_single(dpgmat, dalgmod, nparent, ncross, nprogeny, nhaploblk, objfn_trans, rng):
    opvs = OptimalPopulationValueSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        nhaploblk = nhaploblk,
        objfn_trans = objfn_trans,
        method = "single",
        rng = rng
    )
    miscout = {}
    out_pgmat, out_sel, out_ncross, out_nprogeny = opvs.select(
        pgmat = dpgmat,
        gmat = None,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20,
        miscout = miscout
    )

    assert len(out_sel) == nparent
    assert numpy.all(out_sel < dpgmat.ntaxa)
    assert out_ncross == ncross
    assert out_nprogeny == nprogeny

def test_select_pareto_TypeError(opvs, dpgmat, dalgmod):
    # set to pareto selection
    opvs.method = "pareto"

    # make sure pareto raises errors
    with pytest.raises(TypeError):
        # make selections
        pgmat, sel, ncross, nprogeny = opvs.select(
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
        pgmat, sel, ncross, nprogeny = opvs.select(
            pgmat = dpgmat,
            gmat = None,
            ptdf = None,
            bvmat = None,
            gpmod = None,
            t_cur = 0,
            t_max = 20
        )

def test_select_pareto(nparent, ncross, nprogeny, nhaploblk, method, rng, dpgmat, dalgmod):
    opvs = OptimalPopulationValueSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        nhaploblk = nhaploblk,
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
    sel_pgmat, sel, sel_ncross, sel_nprogeny = opvs.select(
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

def test_pareto(dpgmat, bvmat, dalgmod, nparent, ncross, nprogeny, nhaploblk, rng):
    opvs = OptimalPopulationValueSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        nhaploblk = nhaploblk,
        method = "pareto",
        objfn_trans = None, # sum of two traits
        objfn_trans_kwargs = {}, # no kwargs
        objfn_wt = numpy.array([1.0,1.0]), # maximizing
        ndset_trans = trans_ndpt_to_vec_dist,
        ndset_trans_kwargs = {"objfn_wt": numpy.array([1.0,1.0]), "wt": numpy.array([0.5,0.5])},
        ndset_wt = -1.0,
        rng = rng
    )
    frontier, sel_config = opvs.pareto(
        pgmat = dpgmat,
        gmat = None,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20,
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
    ax.set_title("Optimal Population Value Selection Test Pareto Frontier")
    # ax = pyplot.axes(projection='3d')
    # ax.scatter3D(xdata, ydata, zdata)
    pyplot.savefig("OPV_2d_frontier.png", dpi = 250)

def test_objfn_is_function(dpgmat, bvmat, dalgmod, nparent, ncross, nprogeny, nhaploblk, objfn_trans, rng):
    opvs = OptimalPopulationValueSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        nhaploblk = nhaploblk,
        objfn_trans = objfn_trans,
        method = "single",
        rng = rng
    )
    objfn = opvs.objfn(
        pgmat = dpgmat,
        gmat = None,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )

    assert callable(objfn)

def test_objfn_is_function(opvs, dpgmat, dalgmod):
    # make objective function
    objfn = opvs.objfn(
        pgmat = dpgmat,
        gmat = None,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )
    assert callable(objfn)

def test_objfn_is_function(opvs, dpgmat, dalgmod):
    # make objective function
    objfn = opvs.objfn_vec(
        pgmat = dpgmat,
        gmat = None,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )
    assert callable(objfn)

def test_objfn_is_multiobjective(opvs, dpgmat, dalgmod, mat_ntaxa, mat_ntrait):
    # set transformation function to None
    opvs.objfn_trans = None
    # make objective function
    objfn = opvs.objfn(
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

def test_objfn_vec_is_multiobjective(opvs, dpgmat, dalgmod, mat_ntaxa, mat_ntrait):
    # set transformation function to None
    opvs.objfn_trans = None
    # make objective function
    objfn = opvs.objfn_vec(
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
