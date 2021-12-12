import numpy
import pytest
from numpy.random import Generator
from numpy.random import PCG64
from matplotlib import pyplot

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.breed.prot.sel import ConventionalGenomicSelection
from pybropt.model.gmod import GenericLinearGenomicModel
from pybropt.popgen.gmat import DenseGenotypeMatrix
from pybropt.breed.prot.sel.transfn import trans_ndpt_to_vec_dist
from pybropt.breed.prot.sel.transfn import trans_sum


################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
######################## Genotypes #########################
############################################################
@pytest.fixture
def mat_int8():
    a = numpy.int8([
       [[1, 0, 0, 0, 0, 0, 1, 0, 1, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 1]],
       [[0, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 1, 1, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0]]
    ])
    yield a.sum(0, dtype = 'int8')

@pytest.fixture
def mat_chrgrp():
    yield numpy.int64([1, 1, 2, 2, 3, 3, 4, 4, 5, 5])

@pytest.fixture
def mat_phypos():
    yield numpy.int64([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

@pytest.fixture
def mat_taxa():
    yield numpy.object_(["Line1", "Line2", "Line3", "Line4", "Line5"])

@pytest.fixture
def mat_taxa_grp():
    yield numpy.int64([1, 1, 2, 2, 2])

@pytest.fixture
def dgmat(mat_int8, mat_chrgrp, mat_phypos, mat_taxa, mat_taxa_grp):
    yield DenseGenotypeMatrix(
        mat = mat_int8,
        vrnt_chrgrp = mat_chrgrp,
        vrnt_phypos = mat_phypos,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp
    )

############################################################
###################### Genomic model #######################
############################################################
@pytest.fixture
def beta():
    yield numpy.float64([[1.4, 2.5]])
    # yield numpy.float64([[1.4, 2.5, 7.2]])

@pytest.fixture
def u():
    yield numpy.float64([
        [-0.33,  2.08],
        [-0.69, -1.87],
        [ 1.12,  1.38],
        [-1.44,  0.20],
        [ 0.88, -0.81],
        [ 1.23,  0.25],
        [ 0.19,  4.35],
        [-2.12,  0.73],
        [-0.87,  1.25],
        [ 0.06, -2.52]
    ])
    # yield numpy.float64([
    #     [-0.33,  2.08, -2.42],
    #     [-0.69, -1.87, -1.38],
    #     [ 1.12,  1.38, -5.65],
    #     [-1.44,  0.20,  4.22],
    #     [ 0.88, -0.81,  1.55],
    #     [ 1.23,  0.25,  5.13],
    #     [ 0.19,  4.35,  0.15],
    #     [-2.12,  0.73, -0.38],
    #     [-0.87,  1.25,  2.38],
    #     [ 0.06, -2.52,  2.48]
    # ])

@pytest.fixture
def trait():
    yield numpy.object_(["protein", "yield"])
    # yield numpy.object_(["protein", "yield", "quality"])

@pytest.fixture
def model_name():
    yield "test_glgmod"

@pytest.fixture
def params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def glgmod(beta, u, trait, model_name, params):
    yield GenericLinearGenomicModel(
        beta = beta,
        u = u,
        trait = trait,
        model_name = model_name,
        params = params
    )

############################################################
################## Breeding values model ###################
############################################################
@pytest.fixture
def bvmat(glgmod, dgmat):
    yield glgmod.gebv(dgmat)

############################################################
############### ConventionalGenomicSelection ###############
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
def cgs(nparent, ncross, nprogeny, rng):
    yield ConventionalGenomicSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        rng = rng
    )

@pytest.fixture
def objfn_wt():
    yield [1., 1.]
    # yield [1., 1., 1.]

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(ConventionalGenomicSelection)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(ConventionalGenomicSelection, "__init__")

def test_select_is_concrete():
    generic_assert_concrete_method(ConventionalGenomicSelection, "select")

def test_objfn_is_concrete():
    generic_assert_concrete_method(ConventionalGenomicSelection, "objfn")

def test_objfn_vec_is_concrete():
    generic_assert_concrete_method(ConventionalGenomicSelection, "objfn_vec")

def test_pareto_is_concrete():
    generic_assert_concrete_method(ConventionalGenomicSelection, "pareto")

def test_objfn_static_is_concrete():
    generic_assert_concrete_method(ConventionalGenomicSelection, "objfn_static")

def test_objfn_vec_static_is_concrete():
    generic_assert_concrete_method(ConventionalGenomicSelection, "objfn_vec_static")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_select_single(cgs, dgmat, bvmat, glgmod, nparent, ncross, nprogeny, mat_int8, u):
    pgmat, sel, ncross, nprogeny = cgs.select(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = glgmod,
        t_cur = 0,
        t_max = 20,
        method = "single",
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        objfn_trans = trans_sum,
        objfn_trans_kwargs = {"axis": None},
        objfn_wt = None,
        ndset_trans = None,
        ndset_trans_kwargs = None,
        ndset_wt = None
    )

    assert sel.ndim == 1

def test_select_pareto(cgs, dgmat, bvmat, glgmod, nparent, ncross, nprogeny, objfn_wt):
    pgmat, sel, ncross, nprogeny = cgs.select(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = glgmod,
        t_cur = 0,
        t_max = 20,
        method = "pareto",
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        objfn_trans = None,
        objfn_trans_kwargs = None,
        objfn_wt = objfn_wt,
        ndset_trans = trans_ndpt_to_vec_dist,
        ndset_trans_kwargs = {"objfn_wt": numpy.array(objfn_wt), "wt": numpy.array([.5,.5])},
        ndset_wt = 1.0
    )

    assert sel.ndim == 1

def test_objfn_is_function(cgs, dgmat, bvmat, glgmod):
    objfn = cgs.objfn(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = glgmod,
        t_cur = 0,
        t_max = 20
    )

    assert callable(objfn)

def test_objfn_multiobjective(cgs, dgmat, bvmat, glgmod, mat_int8, u):
    objfn = cgs.objfn(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = glgmod,
        t_cur = 0,
        t_max = 20
    )

    Z = mat_int8        # (n,p) genotypes {0,1,2}
    u = u               # (p,t) regression coefficients
    Y = Z@u             # (n,t) values

    for i,taxon_bv in enumerate(Y):
        numpy.testing.assert_almost_equal(taxon_bv, objfn([i]))

def test_objfn_vec_is_callable(cgs, dgmat, bvmat, glgmod):
    objfn_vec = cgs.objfn_vec(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = glgmod,
        t_cur = 0,
        t_max = 20
    )

    assert callable(objfn_vec)

def test_objfn_vec_multiobjective(cgs, dgmat, bvmat, glgmod, mat_int8, u):
    objfn_vec = cgs.objfn_vec(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = glgmod,
        t_cur = 0,
        t_max = 20
    )

    Z = mat_int8        # (n,p) genotypes {0,1,2}
    u = u               # (p,t) regression coefficients
    Y = Z@u             # (n,t) values

    for i,taxon_bv in enumerate(Y):
        numpy.testing.assert_almost_equal(
            numpy.stack([taxon_bv, taxon_bv]),
            objfn_vec([[i],[i]])
        )

def test_pareto(cgs, dgmat, bvmat, glgmod, objfn_wt):
    frontier, sel_config = cgs.pareto(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = glgmod,
        t_cur = 0,
        t_max = 20,
        objfn_wt = objfn_wt
    )

    xdata = frontier[:,0]
    ydata = frontier[:,1]
    # zdata = frontier[:,2]

    xlabel = glgmod.trait[0]
    ylabel = glgmod.trait[1]

    fig = pyplot.figure()
    ax = pyplot.axes()
    ax.scatter(xdata, ydata)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Conventional Genomic Selection Test Pareto Frontier")
    # ax = pyplot.axes(projection='3d')
    # ax.scatter3D(xdata, ydata, zdata)
    pyplot.savefig("CGS_2d_frontier.png", dpi = 250)

def test_objfn_static_multiobjective(cgs, mat_int8, u):
    Z = mat_int8        # (n,p) genotypes {0,1,2}
    u = u               # (p,t) regression coefficients
    Y = Z@u             # (n,t) values

    for i,taxon_bv in enumerate(Y):
        numpy.testing.assert_almost_equal(
            taxon_bv,
            cgs.objfn_static([i], Z, u, None, None)
        )

def test_objfn_vec_static_multiobjective(cgs, dgmat, bvmat, glgmod, mat_int8, u):
    Z = mat_int8        # (n,p) genotypes {0,1,2}
    u = u               # (p,t) regression coefficients
    Y = Z@u             # (n,t) values

    for i,taxon_bv in enumerate(Y):
        numpy.testing.assert_almost_equal(
            numpy.stack([taxon_bv, taxon_bv]),
            cgs.objfn_vec_static([[i],[i]], Z, u, None, None)
        )
