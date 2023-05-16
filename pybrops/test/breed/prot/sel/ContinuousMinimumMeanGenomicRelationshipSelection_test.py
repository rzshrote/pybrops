import numpy
import pytest
from numpy.random import Generator
from numpy.random import PCG64

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping
from pybrops.breed.prot.sel.transfn import trans_sum
from pybrops.breed.prot.sel.UnconstrainedContinuousMinimumMeanGenomicRelationshipSelection import ContinuousMinimumMeanGenomicRelationshipSelection
from pybrops.popgen.cmat.DenseVanRadenCoancestryMatrix import DenseVanRadenCoancestryMatrix
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
############### ContinuousMinimumMeanGenomicRelationshipSelection ###############
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
def method():
    yield "single"

@pytest.fixture
def ktype():
    yield "gmat"

@pytest.fixture
def kcls():
    yield DenseVanRadenCoancestryMatrix

@pytest.fixture
def mmehs(nparent, ncross, nprogeny, ktype, method, rng):
    yield ContinuousMinimumMeanGenomicRelationshipSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        ktype = ktype,
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
    assert_docstring(ContinuousMinimumMeanGenomicRelationshipSelection)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(ContinuousMinimumMeanGenomicRelationshipSelection, "__init__")

def test_select_is_concrete():
    assert_concrete_method(ContinuousMinimumMeanGenomicRelationshipSelection, "select")

def test_objfn_is_concrete():
    assert_concrete_method(ContinuousMinimumMeanGenomicRelationshipSelection, "objfn")

def test_objfn_vec_is_concrete():
    assert_concrete_method(ContinuousMinimumMeanGenomicRelationshipSelection, "objfn_vec")

def test_pareto_is_concrete():
    assert_concrete_method(ContinuousMinimumMeanGenomicRelationshipSelection, "pareto")

def test_objfn_static_is_concrete():
    assert_concrete_method(ContinuousMinimumMeanGenomicRelationshipSelection, "objfn_static")

def test_objfn_vec_static_is_concrete():
    assert_concrete_method(ContinuousMinimumMeanGenomicRelationshipSelection, "objfn_vec_static")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_select_single(nparent, ncross, nprogeny, ktype, rng, dpgmat, dgmat):
    # make selection object
    selobj = ContinuousMinimumMeanGenomicRelationshipSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        ktype = ktype,
        method = "single",
        objfn_trans = trans_sum, # sum of two traits
        objfn_trans_kwargs = {}, # no kwargs
        objfn_wt = 1.0, # maximizing
        rng = rng
    )

    # set to gmat kinship type
    selobj.ktype = "gmat"
    # make selections
    sel_pgmat, sel, sel_ncross, sel_nprogeny = selobj.select(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    assert id(sel_pgmat) == id(dpgmat)
    assert sel.ndim == 1
    assert sel_ncross == ncross
    assert sel_nprogeny == nprogeny

    # set to pgmat kinship type
    selobj.ktype = "pgmat"
    # make selections
    sel_pgmat, sel, sel_ncross, sel_nprogeny = selobj.select(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    assert id(sel_pgmat) == id(dpgmat)
    assert sel.ndim == 1
    assert sel_ncross == ncross
    assert sel_nprogeny == nprogeny

# # not sure about the purpose of this test
# def test_select_single_objfn_trans_RuntimeError(nparent, ncross, nprogeny, rng, dpgmat, dgmat):
#     # make selection object
#     selobj = ContinuousMinimumMeanGenomicRelationshipSelection(
#         nparent = nparent,
#         ncross = ncross,
#         nprogeny = nprogeny,
#         method = "single",
#         objfn_trans = None, # sum of two traits
#         objfn_trans_kwargs = None, # no kwargs
#         objfn_wt = 1.0, # maximizing function
#         ndset_trans = trans_ndpt_to_vec_dist,
#         ndset_trans_kwargs = {
#             "objfn_wt": numpy.array([1.0, 1.0, 1.0]),   # all objectives maximizing
#             "wt": numpy.array([1./3., 1./3., 1./3.])    # 1/3 equal weight to all
#         },
#         rng = rng
#     )

#     # make sure this raises an error due to lack of shape compatibility
#     with pytest.raises(RuntimeError):
#         # make selections
#         pgmat, sel, ncross, nprogeny = selobj.select(
#             pgmat = dpgmat,
#             gmat = dgmat,
#             ptdf = None,
#             bvmat = None,
#             gpmod = None,
#             t_cur = 0,
#             t_max = 20
#         )

def test_select_pareto_RuntimeError(nparent, ncross, nprogeny, ktype, rng, dpgmat, dgmat):
    # make selection object
    selobj = ContinuousMinimumMeanGenomicRelationshipSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        ktype = ktype,
        method = "pareto",
        objfn_trans = None, # sum of two traits
        objfn_trans_kwargs = None, # no kwargs
        objfn_wt = 1.0, # maximizing function
        rng = rng
    )

    # set to gmat kinship type
    selobj.ktype = "gmat"
    # make sure pareto raises errors
    with pytest.raises(RuntimeError):
        # make selections
        pgmat, sel, ncross, nprogeny = selobj.select(
            pgmat = dpgmat,
            gmat = dgmat,
            ptdf = None,
            bvmat = None,
            gpmod = None,
            t_cur = 0,
            t_max = 20
        )

    # set to gmat kinship type
    selobj.ktype = "pgmat"
    # make sure pareto raises errors
    with pytest.raises(RuntimeError):
        # make selections
        pgmat, sel, ncross, nprogeny = selobj.select(
            pgmat = dpgmat,
            gmat = dgmat,
            ptdf = None,
            bvmat = None,
            gpmod = None,
            t_cur = 0,
            t_max = 20
        )

def test_objfn_is_function(mmehs, dpgmat, dgmat):
    # set to gmat kinship type
    mmehs.ktype = "gmat"
    # make objective function
    objfn = mmehs.objfn(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    assert callable(objfn)

    # set to gmat kinship type
    mmehs.ktype = "pgmat"
    # make objective function
    objfn = mmehs.objfn(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    assert callable(objfn)

def test_objfn_is_multiobjective(nparent, ncross, nprogeny, rng, dpgmat, dgmat):
    # make selection object
    selobj = ContinuousMinimumMeanGenomicRelationshipSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        ktype = "gmat",
        method = "single",
        objfn_trans = None, # sum of two traits
        objfn_trans_kwargs = None, # no kwargs
        objfn_wt = numpy.array([1.0]), # maximizing function
        rng = rng
    )

    # make objective function
    objfn = selobj.objfn(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = None,
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
    assert out.shape[0] == 1

def test_objfn_vec_is_callable(mmehs, dpgmat, dgmat):
    # set to gmat kinship type
    mmehs.ktype = "gmat"
    # make objective function
    objfn_vec = mmehs.objfn_vec(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    assert callable(objfn_vec)

    # set to gmat kinship type
    mmehs.ktype = "pgmat"
    # make objective function
    objfn_vec = mmehs.objfn_vec(
        pgmat = dpgmat,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = None,
        t_cur = 0,
        t_max = 20
    )
    assert callable(objfn_vec)

def test_objfn_vec_is_multiobjective(nparent, ncross, nprogeny, method, rng, dpgmat, dgmat):
    # make selection object
    sel = ContinuousMinimumMeanGenomicRelationshipSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        ktype = "gmat",
        method = "pareto",
        objfn_trans = None, # sum of two traits
        objfn_trans_kwargs = None, # no kwargs
        objfn_wt = numpy.array([1.0]), # maximizing function
        rng = rng
    )

    # make objective function
    objfn = sel.objfn_vec(
        pgmat = None,
        gmat = dgmat,
        ptdf = None,
        bvmat = None,
        gpmod = None,
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
    assert out.shape[1] == 1
