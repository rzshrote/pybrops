import numpy
import pytest
from numpy.random import Generator
from numpy.random import PCG64
from matplotlib import pyplot
from pybrops.breed.prot.sel.EstimatedBreedingValueSelection import EstimatedBreedingValueSubsetSelection
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.breed.prot.sel.prob.trans import trans_ndpt_to_vec_dist
from pybrops.breed.prot.sel.prob.EstimatedBreedingValueSelectionProblem import EstimatedBreedingValueSubsetSelectionProblem

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def ntaxa():
    yield 100

@pytest.fixture
def nvrnt():
    yield 1000

@pytest.fixture
def ntrait():
    yield 2

############################################################
######################## Genotypes #########################
############################################################
@pytest.fixture
def pgmat(ntaxa, nvrnt):
    mat = numpy.random.binomial(1, 0.1, (2,ntaxa,nvrnt)).astype('int8')
    vrnt_chrgrp = numpy.repeat([1,2], nvrnt // 2)
    vrnt_phypos = numpy.arange(1,nvrnt+1)
    vrnt_genpos = numpy.random.random(nvrnt)
    taxa = numpy.array(["Taxa"+str(i).zfill(3) for i in range(1,ntaxa+1)], dtype=object)
    taxa_grp = numpy.repeat([1,2,3,4], ntaxa // 4)
    out = DensePhasedGenotypeMatrix(
        mat = mat,
        vrnt_chrgrp = vrnt_chrgrp,
        vrnt_phypos = vrnt_phypos,
        vrnt_genpos = vrnt_genpos,
        taxa = taxa,
        taxa_grp = taxa_grp
    )
    out.group()
    yield out

@pytest.fixture
def gmat(pgmat):
    dugt = DenseUnphasedGenotyping()
    out = dugt.genotype(pgmat)
    yield out

############################################################
###################### Genomic model #######################
############################################################
@pytest.fixture
def algmod(nvrnt, ntrait):
    beta = numpy.random.uniform(10, 30, (1,ntrait))
    u_misc = None
    u_a = numpy.random.normal(0, 1, (nvrnt,ntrait))
    trait = numpy.array(["Trait"+str(i).zfill(2) for i in range(1,ntrait+1)], dtype=object)
    model_name = "test_dalgmod"
    params = {"a" : 0, "b" : 1}
    yield DenseAdditiveLinearGenomicModel(
        beta = beta,
        u_misc = u_misc,
        u_a = u_a,
        trait = trait,
        model_name = model_name,
        params = params
    )

############################################################
################## Breeding values model ###################
############################################################
@pytest.fixture
def bvmat(algmod, gmat):
    yield algmod.gebv(gmat)

############################################################
######### EstimatedBreedingValueSubsetSelection ##########
############################################################
@pytest.fixture
def nparent():
    yield 5

@pytest.fixture
def ncross():
    yield 1

@pytest.fixture
def nprogeny():
    yield 10

@pytest.fixture
def nobj():
    yield 2

@pytest.fixture
def obj_wt():
    yield numpy.array([1,-1], dtype=float)

@pytest.fixture
def obj_trans():
    yield None

@pytest.fixture
def obj_trans_kwargs():
    yield None

@pytest.fixture
def nineqcv():
    yield 0

@pytest.fixture
def ineqcv_wt():
    yield numpy.array([], dtype=float)

@pytest.fixture
def ineqcv_trans():
    yield None

@pytest.fixture
def ineqcv_trans_kwargs():
    yield None

@pytest.fixture
def neqcv():
    yield 0

@pytest.fixture
def eqcv_wt():
    yield numpy.array([], dtype=float)

@pytest.fixture
def eqcv_trans():
    yield None

@pytest.fixture
def eqcv_trans_kwargs():
    yield None

@pytest.fixture
def ndset_wt():
    yield -1.0

@pytest.fixture
def ndset_trans():
    yield trans_ndpt_to_vec_dist

@pytest.fixture
def ndset_trans_kwargs(obj_wt):
    yield {"obj_wt": obj_wt, "vec_wt": numpy.array([0.5,0.5])}

@pytest.fixture
def method():
    yield "single"

@pytest.fixture
def rng():
    yield Generator(PCG64(192837465))

@pytest.fixture
def soalgo():
    yield None

@pytest.fixture
def moalgo():
    yield None

@pytest.fixture
def selprot(
    nparent, ncross, nprogeny, method,
    nobj, obj_wt, obj_trans, obj_trans_kwargs, 
    nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs, 
    neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs, 
    ndset_wt, ndset_trans,  ndset_trans_kwargs, 
    rng, soalgo, moalgo
):
    yield EstimatedBreedingValueSubsetSelection(
        nparent = 5, 
        nmating = 1, 
        nprogeny = 10,
        method = method,
        nobj = nobj,
        obj_wt = obj_wt,
        obj_trans = obj_trans,
        obj_trans_kwargs = obj_trans_kwargs,
        nineqcv = nineqcv,
        ineqcv_wt = ineqcv_wt,
        ineqcv_trans = ineqcv_trans,
        ineqcv_trans_kwargs = ineqcv_trans_kwargs,
        neqcv = neqcv,
        eqcv_wt = eqcv_wt,
        eqcv_trans = eqcv_trans,
        eqcv_trans_kwargs = eqcv_trans_kwargs,
        ndset_wt = ndset_wt,
        ndset_trans = ndset_trans, 
        ndset_trans_kwargs = ndset_trans_kwargs, 
        rng = rng, 
        soalgo = soalgo,
        moalgo = moalgo
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(EstimatedBreedingValueSubsetSelection)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(EstimatedBreedingValueSubsetSelection, "__init__")

def test_objfn_is_concrete():
    assert_concrete_method(EstimatedBreedingValueSubsetSelection, "problem")

def test_pareto_is_concrete():
    assert_concrete_method(EstimatedBreedingValueSubsetSelection, "pareto")

def test_select_is_concrete():
    assert_concrete_method(EstimatedBreedingValueSubsetSelection, "select")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################

def test_problem(selprot, pgmat, gmat, bvmat, algmod):
    prob = selprot.problem(
        pgmat = pgmat,
        gmat = gmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = algmod,
        t_cur = None,
        t_max = None
    )

    assert isinstance(prob, SelectionProblem)
    assert isinstance(prob, SubsetSelectionProblem)
    assert isinstance(prob, EstimatedBreedingValueSubsetSelectionProblem)


def test_pareto(selprot, pgmat, gmat, bvmat, algmod):
    frontier, sel_config = selprot.pareto(
        pgmat = pgmat,
        gmat = gmat,
        ptdf = None,
        bvmat = bvmat,
        gpmod = algmod,
        t_cur = 0,
        t_max = 20
    )

    xdata = frontier[:,0]
    ydata = frontier[:,1]
    # zdata = frontier[:,2]

    xlabel = algmod.trait[0]
    ylabel = algmod.trait[1]

    fig = pyplot.figure()
    ax = pyplot.axes()
    ax.scatter(xdata, ydata)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("EBV Subset Selection Test Pareto Frontier")
    # ax = pyplot.axes(projection='3d')
    # ax.scatter3D(xdata, ydata, zdata)
    pyplot.savefig("EBV_Subset_2d_frontier.png", dpi = 250)

    print(sel_config)
    print(sel_config.sum(0))
    print(sel_config.sum(1))
    assert False

# def test_select_single(selprot, pgmat, gmat, bvmat, algmod):
#     selprot.method = "single"

#     pgmat, sel, ncross, nprogeny = selprot.select(
#         pgmat = pgmat,
#         gmat = gmat,
#         ptdf = None,
#         bvmat = bvmat,
#         gpmod = algmod,
#         t_cur = 0,
#         t_max = 20,
#         miscout = None
#     )

#     assert sel.ndim == 1

# def test_select_pareto(selprot, pgmat, gmat, bvmat, algmod):
#     selprot.method = "pareto"

#     pgmat, sel, ncross, nprogeny = selprot.select(
#         pgmat = pgmat,
#         gmat = gmat,
#         ptdf = None,
#         bvmat = bvmat,
#         gpmod = algmod,
#         t_cur = 0,
#         t_max = 20,
#         miscout = None
#     )

#     assert sel.ndim == 1

