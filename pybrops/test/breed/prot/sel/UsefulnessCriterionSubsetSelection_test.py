import pytest
from pybrops.breed.prot.sel.UsefulnessCriterionSelection import UsefulnessCriterionSubsetSelection
from pybrops.breed.prot.sel.prob.SubsetMateSelectionProblem import SubsetMateSelectionProblem
from pybrops.test.breed.prot.sel.common_fixtures_large import *
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method


################################ Test fixtures #################################

@pytest.fixture
def selprot(
        common_ntrait,
        common_nself,
        common_upper_percentile,
        common_vmatfcty,
        common_gmapfn,
        common_unique_parents,
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_nobj,
        common_obj_wt,
        common_obj_trans,
        common_obj_trans_kwargs,
        common_nineqcv,
        common_ineqcv_wt,
        common_ineqcv_trans,
        common_ineqcv_trans_kwargs,
        common_neqcv,
        common_eqcv_wt,
        common_eqcv_trans,
        common_eqcv_trans_kwargs,
        common_ndset_wt,
        common_ndset_trans,
        common_ndset_trans_kwargs,
        common_rng,
        common_soalgo,
        common_moalgo
    ):
    out = UsefulnessCriterionSubsetSelection(
        ntrait = common_ntrait,
        nself = common_nself,
        upper_percentile = common_upper_percentile,
        vmatfcty = common_vmatfcty,
        gmapfn = common_gmapfn,
        unique_parents = common_unique_parents,
        ncross = common_ncross,
        nparent = common_nparent,
        nmating = common_nmating,
        nprogeny = common_nprogeny,
        nobj = common_nobj,
        obj_wt = common_obj_wt,
        obj_trans = common_obj_trans,
        obj_trans_kwargs = common_obj_trans_kwargs,
        nineqcv = common_nineqcv,
        ineqcv_wt = common_ineqcv_wt,
        ineqcv_trans = common_ineqcv_trans,
        ineqcv_trans_kwargs = common_ineqcv_trans_kwargs,
        neqcv = common_neqcv,
        eqcv_wt = common_eqcv_wt,
        eqcv_trans = common_eqcv_trans,
        eqcv_trans_kwargs = common_eqcv_trans_kwargs,
        ndset_wt = common_ndset_wt,
        ndset_trans = common_ndset_trans,
        ndset_trans_kwargs = common_ndset_trans_kwargs,
        rng = common_rng,
        soalgo = common_soalgo,
        moalgo = common_moalgo
    )
    yield out

############################## Test class docstring ############################
def test_class_docstring():
    assert_docstring(UsefulnessCriterionSubsetSelection)

############################# Test concrete methods ############################
def test_init_is_concrete():
    assert_concrete_method(UsefulnessCriterionSubsetSelection, "__init__")

def test_problem_is_concrete():
    assert_concrete_method(UsefulnessCriterionSubsetSelection, "problem")

def test_sosolve_is_concrete():
    assert_concrete_method(UsefulnessCriterionSubsetSelection, "sosolve")

def test_mosolve_is_concrete():
    assert_concrete_method(UsefulnessCriterionSubsetSelection, "mosolve")

def test_select_is_concrete():
    assert_concrete_method(UsefulnessCriterionSubsetSelection, "select")

###################### Test concrete method functionality ######################
def test_problem(
        selprot,
        common_pgmat,
        common_gmat,
        common_ptdf,
        common_bvmat,
        common_gpmod,
        common_t_cur,
        common_t_max
    ):
    # create selection problem
    prob = selprot.problem(
        pgmat = common_pgmat,
        gmat = common_gmat,
        ptdf = common_ptdf,
        bvmat = common_bvmat,
        gpmod = common_gpmod,
        t_cur = common_t_cur,
        t_max = common_t_max
    )

    # check that it is the right type
    assert isinstance(prob, SubsetMateSelectionProblem)