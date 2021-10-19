import numpy
import pytest
from numpy.random import Generator
from numpy.random import PCG64

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
    yield numpy.float64([[1.4, 2.5, 7.2]])

@pytest.fixture
def u():
    yield numpy.float64([
        [-0.33,  2.08, -2.42],
        [-0.69, -1.87, -1.38],
        [ 1.12,  1.38, -5.65],
        [-1.44,  0.20,  4.22],
        [ 0.88, -0.81,  1.55],
        [ 1.23,  0.25,  5.13],
        [ 0.19,  4.35,  0.15],
        [-2.12,  0.73, -0.38],
        [-0.87,  1.25,  2.38],
        [ 0.06, -2.52,  2.48]
    ])

@pytest.fixture
def trait():
    yield numpy.object_(["protein", "yield", "quality"])

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

# TODO:
# def test_pareto_is_concrete():
#     generic_assert_concrete_method(ConventionalGenomicSelection, "pareto")

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
def test_objfn_multiobjective(cgs, dgmat, bvmat, glgmod, ncross, nprogeny, mat_int8, u):
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

    Z = mat_int8        # (n,p) genotypes {0,1,2}
    u = u               # (p,t) regression coefficients
    Y = Z@u             # (n,t) values

    for i,taxon_bv in enumerate(Y):
        numpy.testing.assert_almost_equal(taxon_bv, objfn([i]))
        # print(taxon_bv == objfn([i]))
        # assert numpy.all(taxon_bv == objfn([i]))

# def test_pselect(cgs, dgmat, bvmat, glgmod, ncross, nprogeny):
#     a,b,c,d,e = cgs.select(
#         pgmat = None,
#         gmat = dgmat,
#         ptdf = None,
#         bvmat = bvmat,
#         gpmod = glgmod,
#         t_cur = 0,
#         t_max = 20
#     )
#
#     assert numpy.all(b == [3,4])
#     assert c == ncross
#     assert d == nprogeny
