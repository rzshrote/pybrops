import numpy
import pytest

from numpy.random import Generator
from numpy.random import PCG64
from matplotlib import pyplot

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.prot.sel.UnconstrainedMultiObjectiveGenomicSelection import MultiObjectiveGenomicSelection
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.bvmat.DenseEstimatedBreedingValueMatrix import DenseEstimatedBreedingValueMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.opt.algo.UnconstrainedSteepestAscentSetHillClimber import UnconstrainedSteepestAscentSetHillClimber
from pybrops.opt.algo.UnconstrainedStochasticAscentSetHillClimber import UnconstrainedStochasticAscentSetHillClimber
from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.breed.prot.sel.transfn import trans_sum
from pybrops.breed.prot.sel.transfn import trans_ndpt_to_vec_dist

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
def beta():
    yield numpy.float64([
        [25.6]
    ])

@pytest.fixture
def mat_u_misc():
    yield None

@pytest.fixture
def mat_u_a():
    yield numpy.float64([
        [ 0.17],
        [-0.51],
        [ 1.99],
        [-0.33],
        [ 0.28],
        [-2.16],
        [-1.  ],
        [-0.99],
        [ 0.76],
        [-2.37],
        [ 0.56],
        [ 0.2 ],
        [-0.  ],
        [-1.21],
        [-1.24],
        [-0.21],
        [-0.84],
        [-1.3 ],
        [-1.09],
        [ 0.71]
    ])

@pytest.fixture
def trait():
    yield numpy.object_(["yield"])

@pytest.fixture
def model_name():
    yield "test_glgmod"

@pytest.fixture
def params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def dalgmod(beta, mat_u_misc, mat_u_a, trait, model_name, params):
    yield DenseAdditiveLinearGenomicModel(
        beta = beta,
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
############# MultiObjectiveGenomicSelection ##############
############################################################
@pytest.fixture
def nparent():
    yield 6

@pytest.fixture
def ncross():
    yield 1

@pytest.fixture
def nprogeny():
    yield 10

@pytest.fixture
def objfn_trans():
    yield trans_sum

@pytest.fixture
def objfn_wt():
    yield -1.0 # since this a minimizing function

@pytest.fixture
def rng():
    yield Generator(PCG64(2367929260))

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(MultiObjectiveGenomicSelection)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(MultiObjectiveGenomicSelection, "__init__")

def test_select_is_concrete():
    assert_concrete_method(MultiObjectiveGenomicSelection, "select")

def test_objfn_is_concrete():
    assert_concrete_method(MultiObjectiveGenomicSelection, "objfn")

def test_objfn_vec_is_concrete():
    assert_concrete_method(MultiObjectiveGenomicSelection, "objfn_vec")

def test_pareto_is_concrete():
    assert_concrete_method(MultiObjectiveGenomicSelection, "pareto")

def test_objfn_static_is_concrete():
    assert_concrete_method(MultiObjectiveGenomicSelection, "objfn_static")

def test_objfn_vec_static_is_concrete():
    assert_concrete_method(MultiObjectiveGenomicSelection, "objfn_vec_static")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################

########################################
########### test calc_mkrwt ############
########################################
def test_calc_mkrwt_magnitude(mat_u_a):
    a = MultiObjectiveGenomicSelection._calc_mkrwt("magnitude", mat_u_a)
    b = numpy.absolute(mat_u_a)
    assert numpy.all(a == b)

def test_calc_mkrwt_equal(mat_u_a):
    a = MultiObjectiveGenomicSelection._calc_mkrwt("equal", mat_u_a)
    assert numpy.all(a == 1.0)

def test_calc_mkrwt_str_case(mat_u_a):
    a = MultiObjectiveGenomicSelection._calc_mkrwt("mAgNiTuDe", mat_u_a)
    b = numpy.absolute(mat_u_a)
    assert numpy.all(a == b)
    a = MultiObjectiveGenomicSelection._calc_mkrwt("Equal", mat_u_a)
    assert numpy.all(a == 1.0)

def test_calc_mkrwt_str_ValueError(mat_u_a):
    with pytest.raises(ValueError):
        a = MultiObjectiveGenomicSelection._calc_mkrwt("unknown", mat_u_a)

def test_calc_mkrwt_ndarray(mat_u_a):
    wt = numpy.random.normal(size = mat_u_a.shape)
    a = MultiObjectiveGenomicSelection._calc_mkrwt(wt, mat_u_a)
    assert numpy.all(a == wt)

def test_calc_mkrwt_type_TypeError(mat_u_a):
    with pytest.raises(TypeError):
        a = MultiObjectiveGenomicSelection._calc_mkrwt(None, mat_u_a)

########################################
########### test calc_tfreq ############
########################################
def test_calc_tfreq_positive(mat_u_a):
    a = MultiObjectiveGenomicSelection._calc_tfreq("positive", mat_u_a)
    b = numpy.float64(mat_u_a >= 0.0)
    assert numpy.all(a == b)

def test_calc_tfreq_negative(mat_u_a):
    a = MultiObjectiveGenomicSelection._calc_tfreq("negative", mat_u_a)
    b = numpy.float64(mat_u_a <= 0.0)
    assert numpy.all(a == b)

def test_calc_tfreq_stabilizing(mat_u_a):
    a = MultiObjectiveGenomicSelection._calc_tfreq("stabilizing", mat_u_a)
    assert numpy.all(a == 0.5)

def test_calc_tfreq_str_case(mat_u_a):
    a = MultiObjectiveGenomicSelection._calc_tfreq("PoSiTiVe", mat_u_a)
    b = numpy.float64(mat_u_a >= 0.0)
    assert numpy.all(a == b)
    a = MultiObjectiveGenomicSelection._calc_tfreq("NEGATIVE", mat_u_a)
    b = numpy.float64(mat_u_a <= 0.0)
    assert numpy.all(a == b)
    a = MultiObjectiveGenomicSelection._calc_tfreq("Stabilizing", mat_u_a)
    assert numpy.all(a == 0.5)

def test_calc_tfreq_str_ValueError(mat_u_a):
    with pytest.raises(ValueError):
        a = MultiObjectiveGenomicSelection._calc_tfreq("unknown", mat_u_a)

def test_calc_tfreq_ndarray(mat_u_a):
    wt = numpy.random.uniform(0, 1, size = mat_u_a.shape)
    a = MultiObjectiveGenomicSelection._calc_tfreq(wt, mat_u_a)
    assert numpy.all(a == wt)

def test_calc_tfreq_type_TypeError(mat_u_a):
    with pytest.raises(TypeError):
        a = MultiObjectiveGenomicSelection._calc_tfreq(None, mat_u_a)

########################################
######### selection functions ##########
########################################
def test_select_single(dpgmat, dgmat, bvmat, dalgmod, nparent, ncross, nprogeny, objfn_trans, rng):
    mogs = MultiObjectiveGenomicSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        objfn_trans = objfn_trans,
        method = "single",
        rng = rng
    )
    miscout = {}
    out_pgmat, out_sel, out_ncross, out_nprogeny = mogs.select(
        pgmat = dpgmat,
        gmat = dgmat,
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

def test_objfn_is_function(dpgmat, dgmat, bvmat, dalgmod, nparent, ncross, nprogeny, objfn_trans, rng):
    mogs = MultiObjectiveGenomicSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        unique_parents = True,
        objfn_trans = objfn_trans,
        method = "single",
        rng = rng
    )
    objfn = mogs.objfn(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )

    assert callable(objfn)

def test_pareto(dpgmat, dgmat, bvmat, dalgmod, nparent, ncross, nprogeny, rng):
    mogs = MultiObjectiveGenomicSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        unique_parents = True,
        objfn_trans = None,
        objfn_trans_kwargs = None,
        objfn_wt = numpy.array([-1.0, -1.0]),
        rng = rng
    )
    frontier, sel_config = mogs.pareto(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20,
    )

    xdata = frontier[:,0]
    ydata = frontier[:,1]
    # zdata = frontier[:,2]

    xlabel = dalgmod.trait[0] + " PAU"
    ylabel = dalgmod.trait[0] + " PAFD"

    fig = pyplot.figure()
    ax = pyplot.axes()
    ax.scatter(xdata, ydata)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Multi-Objective Genomic Selection Test Pareto Frontier")
    # ax = pyplot.axes(projection='3d')
    # ax.scatter3D(xdata, ydata, zdata)
    pyplot.savefig("MOGS_2d_frontier.png", dpi = 250)
