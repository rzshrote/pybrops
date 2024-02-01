import os
import pytest
from matplotlib import pyplot
from pybrops.breed.prot.sel.ExpectedMaximumBreedingValueSelection import ExpectedMaximumBreedingValueSubsetSelection
from pybrops.breed.prot.sel.cfg.SubsetSelectionConfiguration import SubsetSelectionConfiguration
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pybrops.breed.prot.sel.soln.SubsetSelectionSolution import SubsetSelectionSolution
from .common_fixtures_large import *
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method

################################ Test fixtures #################################

@pytest.fixture
def selprot_single(
        common_ntrait,
        common_nrep,
        common_mateprot,
        common_unique_parents,
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_nobj_single,
        common_obj_wt_single,
        common_obj_trans_single,
        common_obj_trans_kwargs_single,
        common_nineqcv_single,
        common_ineqcv_wt_single,
        common_ineqcv_trans_single,
        common_ineqcv_trans_kwargs_single,
        common_neqcv_single,
        common_eqcv_wt_single,
        common_eqcv_trans_single,
        common_eqcv_trans_kwargs_single,
        common_ndset_wt_single,
        common_ndset_trans_single,
        common_ndset_trans_kwargs_single,
        common_rng,
        common_soalgo,
        common_moalgo
    ):
    out = ExpectedMaximumBreedingValueSubsetSelection(
        ntrait = common_ntrait,
        nrep = common_nrep,
        mateprot = common_mateprot,
        unique_parents = common_unique_parents,
        ncross = common_ncross,
        nparent = common_nparent,
        nmating = common_nmating,
        nprogeny = common_nprogeny,
        nobj = common_nobj_single,
        obj_wt = common_obj_wt_single,
        obj_trans = common_obj_trans_single,
        obj_trans_kwargs = common_obj_trans_kwargs_single,
        nineqcv = common_nineqcv_single,
        ineqcv_wt = common_ineqcv_wt_single,
        ineqcv_trans = common_ineqcv_trans_single,
        ineqcv_trans_kwargs = common_ineqcv_trans_kwargs_single,
        neqcv = common_neqcv_single,
        eqcv_wt = common_eqcv_wt_single,
        eqcv_trans = common_eqcv_trans_single,
        eqcv_trans_kwargs = common_eqcv_trans_kwargs_single,
        ndset_wt = common_ndset_wt_single,
        ndset_trans = common_ndset_trans_single,
        ndset_trans_kwargs = common_ndset_trans_kwargs_single,
        rng = common_rng,
        soalgo = common_soalgo,
        moalgo = common_moalgo
    )
    yield out

@pytest.fixture
def selprot_multi(
        common_ntrait,
        common_nrep,
        common_mateprot,
        common_unique_parents,
        common_ncross,
        common_nparent,
        common_nmating,
        common_nprogeny,
        common_nobj_multi,
        common_obj_wt_multi,
        common_obj_trans_multi,
        common_obj_trans_kwargs_multi,
        common_nineqcv_multi,
        common_ineqcv_wt_multi,
        common_ineqcv_trans_multi,
        common_ineqcv_trans_kwargs_multi,
        common_neqcv_multi,
        common_eqcv_wt_multi,
        common_eqcv_trans_multi,
        common_eqcv_trans_kwargs_multi,
        common_ndset_wt_multi,
        common_ndset_trans_multi,
        common_ndset_trans_kwargs_multi,
        common_rng,
        common_soalgo,
        common_moalgo
    ):
    out = ExpectedMaximumBreedingValueSubsetSelection(
        ntrait = common_ntrait,
        nrep = common_nrep,
        mateprot = common_mateprot,
        unique_parents = common_unique_parents,
        ncross = common_ncross,
        nparent = common_nparent,
        nmating = common_nmating,
        nprogeny = common_nprogeny,
        nobj = common_nobj_multi,
        obj_wt = common_obj_wt_multi,
        obj_trans = common_obj_trans_multi,
        obj_trans_kwargs = common_obj_trans_kwargs_multi,
        nineqcv = common_nineqcv_multi,
        ineqcv_wt = common_ineqcv_wt_multi,
        ineqcv_trans = common_ineqcv_trans_multi,
        ineqcv_trans_kwargs = common_ineqcv_trans_kwargs_multi,
        neqcv = common_neqcv_multi,
        eqcv_wt = common_eqcv_wt_multi,
        eqcv_trans = common_eqcv_trans_multi,
        eqcv_trans_kwargs = common_eqcv_trans_kwargs_multi,
        ndset_wt = common_ndset_wt_multi,
        ndset_trans = common_ndset_trans_multi,
        ndset_trans_kwargs = common_ndset_trans_kwargs_multi,
        rng = common_rng,
        soalgo = common_soalgo,
        moalgo = common_moalgo
    )
    yield out

############################## Test class docstring ############################
def test_class_docstring():
    assert_docstring(ExpectedMaximumBreedingValueSubsetSelection)

############################# Test concrete methods ############################
def test_init_is_concrete():
    assert_concrete_method(ExpectedMaximumBreedingValueSubsetSelection, "__init__")

def test_problem_is_concrete():
    assert_concrete_method(ExpectedMaximumBreedingValueSubsetSelection, "problem")

def test_sosolve_is_concrete():
    assert_concrete_method(ExpectedMaximumBreedingValueSubsetSelection, "sosolve")

def test_mosolve_is_concrete():
    assert_concrete_method(ExpectedMaximumBreedingValueSubsetSelection, "mosolve")

def test_select_is_concrete():
    assert_concrete_method(ExpectedMaximumBreedingValueSubsetSelection, "select")

###################### Test concrete method functionality ######################
def test_problem(
        selprot_single,
        selprot_multi,
        common_pgmat,
        common_gmat,
        common_ptdf,
        common_bvmat,
        common_gpmod,
        common_t_cur,
        common_t_max
    ):
    # create selection problem
    prob = selprot_single.problem(
        pgmat = common_pgmat,
        gmat = common_gmat,
        ptdf = common_ptdf,
        bvmat = common_bvmat,
        gpmod = common_gpmod,
        t_cur = common_t_cur,
        t_max = common_t_max
    )

    # check that it is the right type
    assert isinstance(prob, SubsetSelectionProblem)

    # create selection problem
    prob = selprot_multi.problem(
        pgmat = common_pgmat,
        gmat = common_gmat,
        ptdf = common_ptdf,
        bvmat = common_bvmat,
        gpmod = common_gpmod,
        t_cur = common_t_cur,
        t_max = common_t_max
    )

    # check that it is the right type
    assert isinstance(prob, SubsetSelectionProblem)

def test_sosolve(
        selprot_single,
        selprot_multi,
        common_pgmat,
        common_gmat,
        common_ptdf,
        common_bvmat,
        common_gpmod,
        common_t_cur,
        common_t_max
    ):
    # solve single objective problem
    soln = selprot_single.sosolve(
        common_pgmat,
        common_gmat,
        common_ptdf,
        common_bvmat,
        common_gpmod,
        common_t_cur,
        common_t_max
    )

    # test for right type
    assert isinstance(soln, SubsetSelectionSolution)

    # make sure multi objective problem raises error
    with pytest.raises(RuntimeError):
        soln = selprot_multi.sosolve(
            common_pgmat,
            common_gmat,
            common_ptdf,
            common_bvmat,
            common_gpmod,
            common_t_cur,
            common_t_max
        )

def test_mosolve(
        selprot_single,
        selprot_multi,
        common_pgmat,
        common_gmat,
        common_ptdf,
        common_bvmat,
        common_gpmod,
        common_t_cur,
        common_t_max
    ):
    # solve single objective problem
    soln = selprot_multi.mosolve(
        common_pgmat,
        common_gmat,
        common_ptdf,
        common_bvmat,
        common_gpmod,
        common_t_cur,
        common_t_max
    )

    # test for right type
    assert isinstance(soln, SubsetSelectionSolution)

    # make a directory if needed
    if not os.path.isdir("frontier_plots"):
        os.makedirs("frontier_plots")
    
    # plot data
    fig = pyplot.figure()
    ax = pyplot.axes()
    ax.scatter(
        soln.soln_obj[:,0], 
        soln.soln_obj[:,1]
    )
    ax.set_xlabel("Trait 1")
    ax.set_ylabel("Trait 2")
    ax.set_title("Random Subset Selection Test Pareto Frontier")
    pyplot.savefig("frontier_plots/ExpectedMaximumBreedingValueSubsetSelection_2d_frontier.png", dpi = 250)
    pyplot.close()

    # make sure multi objective problem raises error
    with pytest.raises(RuntimeError):
        soln = selprot_single.mosolve(
            common_pgmat,
            common_gmat,
            common_ptdf,
            common_bvmat,
            common_gpmod,
            common_t_cur,
            common_t_max
        )

def test_select(
        selprot_single,
        selprot_multi,
        common_pgmat,
        common_gmat,
        common_ptdf,
        common_bvmat,
        common_gpmod,
        common_t_cur,
        common_t_max
    ):
    # make single-objective selections
    selcfg_single = selprot_single.select(
        common_pgmat,
        common_gmat,
        common_ptdf,
        common_bvmat,
        common_gpmod,
        common_t_cur,
        common_t_max
    )

    assert isinstance(selcfg_single, SubsetSelectionConfiguration)

    # make multi-objective selections
    selcfg_multi = selprot_multi.select(
        common_pgmat,
        common_gmat,
        common_ptdf,
        common_bvmat,
        common_gpmod,
        common_t_cur,
        common_t_max
    )

    assert isinstance(selcfg_multi, SubsetSelectionConfiguration)
