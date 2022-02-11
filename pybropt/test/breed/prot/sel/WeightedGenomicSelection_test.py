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

from pybropt.breed.prot.sel.WeightedGenomicSelection import WeightedGenomicSelection
from pybropt.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybropt.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybropt.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybropt.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping


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
############### WeightedGenomicSelection ###############
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
def wgs(nparent, ncross, nprogeny, rng):
    yield WeightedGenomicSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        rng = rng
    )

@pytest.fixture
def objfn_wt():
    yield [1., 1.]

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(WeightedGenomicSelection)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(WeightedGenomicSelection, "__init__")

def test_select_is_concrete():
    generic_assert_concrete_method(WeightedGenomicSelection, "select")

def test_objfn_is_concrete():
    generic_assert_concrete_method(WeightedGenomicSelection, "objfn")

def test_objfn_vec_is_concrete():
    generic_assert_concrete_method(WeightedGenomicSelection, "objfn_vec")

def test_pareto_is_concrete():
    generic_assert_concrete_method(WeightedGenomicSelection, "pareto")

def test_objfn_static_is_concrete():
    generic_assert_concrete_method(WeightedGenomicSelection, "objfn_static")

def test_objfn_vec_static_is_concrete():
    generic_assert_concrete_method(WeightedGenomicSelection, "objfn_vec_static")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_objfn_multiobjective(wgs, dgmat, bvmat, dalgmod, ncross, nprogeny, mat_int8, mat_u_a):
    objfn = wgs.objfn(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20
    )

    assert callable(objfn)

    geno = mat_int8.sum(0)
    afreq = ((1/(geno.shape[0]*2)) * geno.sum(0))[:,None]
    fafreq = numpy.where(               # (p,t) calculate favorable allele frequencies
        mat_u_a > 0.0,                        # if dominant (1) allele is beneficial
        afreq,                          # get dominant allele frequency
        1.0 - afreq                     # else get recessive allele frequency
    )
    fafreq[fafreq <= 0.0] = 1.0         # avoid division by zero/imaginary
    uwt = numpy.power(fafreq, -0.5)  # calculate weights: 1/sqrt(p)

    Z = geno            # (n,p) genotypes {0,1,2}
    u_a = mat_u_a       # (p,t) regression coefficients
    Y = Z@(u_a*uwt)     # (n,t) values

    for i,taxon_bv in enumerate(Y):
        numpy.testing.assert_almost_equal(taxon_bv, objfn([i]))
        # print(taxon_bv == objfn([i]))
        # assert numpy.all(taxon_bv == objfn([i]))

def test_pareto(wgs, dgmat, bvmat, dalgmod, objfn_wt):
    frontier, sel_config = wgs.pareto(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20,
        objfn_wt = objfn_wt
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
    ax.set_title("Weighted Genomic Selection Test Pareto Frontier")
    # ax = pyplot.axes(projection='3d')
    # ax.scatter3D(xdata, ydata, zdata)
    pyplot.savefig("WGS_2d_frontier.png", dpi = 250)


# def test_pselect(wgs, dgmat, bvmat, dalgmod, ncross, nprogeny):
#     a,b,c,d,e = wgs.select(
#         pgmat = None,
#         gmat = dgmat,
#         ptdf = None,
#         bvmat = bvmat,
#         gpmod = dalgmod,
#         t_cur = 0,
#         t_max = 20
#     )
#
#     assert numpy.all(b == [3,4])
#     assert c == ncross
#     assert d == nprogeny
