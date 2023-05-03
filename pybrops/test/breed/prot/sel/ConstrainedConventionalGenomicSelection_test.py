import numpy
import pytest
from numpy.random import Generator
from numpy.random import PCG64
from matplotlib import pyplot
from pybrops.breed.prot.sel.prob.SelectionProblemType import SelectionProblemType

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_abstract_function
from pybrops.test import assert_abstract_property
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.breed.prot.sel.ConstrainedConventionalGenomicSelection import ConstrainedConventionalGenomicSelection
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.breed.prot.sel.prob.trans import trans_ndpt_to_vec_dist

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
def mat_int8(ntaxa, nvrnt):
    yield numpy.random.binomial(1, 0.1, (2,ntaxa,nvrnt)).astype('int8')

@pytest.fixture
def mat_chrgrp(nvrnt):
    yield numpy.repeat([1,2], nvrnt // 2)

@pytest.fixture
def mat_phypos(nvrnt):
    yield numpy.arange(1,nvrnt+1)

@pytest.fixture
def mat_genpos(nvrnt):
    out = numpy.random.random(nvrnt)
    out.sort()
    yield out

@pytest.fixture
def mat_taxa(ntaxa):
    yield numpy.array(["Line"+str(i).zfill(3) for i in range(1,ntaxa+1)], dtype=object)

@pytest.fixture
def mat_taxa_grp(ntaxa):
    yield numpy.repeat([1,2,3,4], ntaxa // 4)

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
def mat_beta(ntrait):
    yield numpy.random.uniform(10, 30, (1,ntrait))

@pytest.fixture
def mat_u_misc():
    yield None

@pytest.fixture
def mat_u_a(nvrnt, ntrait):
    yield numpy.random.normal(0, 1, (nvrnt,ntrait))

@pytest.fixture
def trait(ntrait):
    yield numpy.array(["Trait"+str(i).zfill(2) for i in range(1,ntrait+1)], dtype=object)

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
######### ConstrainedConventionalGenomicSelection ##########
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
def cgs(
    nparent, 
    ncross, 
    nprogeny,
    nobj,
    obj_wt,
    obj_trans,
    obj_trans_kwargs,
    nineqcv,
    ineqcv_wt,
    ineqcv_trans,
    ineqcv_trans_kwargs,
    neqcv,
    eqcv_wt,
    eqcv_trans,
    eqcv_trans_kwargs,
    ndset_wt,
    ndset_trans, 
    ndset_trans_kwargs, 
    method,
    rng, 
    soalgo,
    moalgo
):
    yield ConstrainedConventionalGenomicSelection(
        nparent = nparent, 
        ncross = ncross, 
        nprogeny = nprogeny,
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
        method = method,
        rng = rng, 
        soalgo = soalgo,
        moalgo = moalgo
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(ConstrainedConventionalGenomicSelection)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(ConstrainedConventionalGenomicSelection, "__init__")

def test_objfn_is_concrete():
    assert_concrete_method(ConstrainedConventionalGenomicSelection, "problem")

def test_pareto_is_concrete():
    assert_concrete_method(ConstrainedConventionalGenomicSelection, "pareto")

def test_select_is_concrete():
    assert_concrete_method(ConstrainedConventionalGenomicSelection, "select")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################

def test_problem(cgs, dgmat, dalgmod):
    prob = cgs.problem(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = None,
        t_max = None
    )

    assert isinstance(prob, SelectionProblemType)

def test_pareto(cgs, dgmat, dalgmod):
    frontier, sel_config = cgs.pareto(
        pgmat = None,
        gmat = dgmat,
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
    ax.set_title("Constrained Conventional Genomic Selection Test Pareto Frontier")
    # ax = pyplot.axes(projection='3d')
    # ax.scatter3D(xdata, ydata, zdata)
    pyplot.savefig("CGS_2d_frontier.png", dpi = 250)

def test_select_single(cgs, dgmat, dalgmod):
    cgs.method = "single"

    pgmat, sel, ncross, nprogeny = cgs.select(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20,
        miscout = None
    )

    assert sel.ndim == 1

def test_select_pareto(cgs, dgmat, dalgmod):
    cgs.method = "pareto"

    pgmat, sel, ncross, nprogeny = cgs.select(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = dalgmod,
        t_cur = 0,
        t_max = 20,
        miscout = None
    )

    assert sel.ndim == 1

