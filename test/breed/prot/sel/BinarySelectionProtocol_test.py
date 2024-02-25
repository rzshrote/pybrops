import pytest
from pybrops.breed.prot.sel.BinarySelectionProtocol import BinarySelectionProtocol
from pybrops.breed.prot.sel.BinarySelectionProtocol import check_is_BinarySelectionProtocol
from pybrops.test.assert_python import assert_method_isabstract, assert_method_isconcrete, assert_class_documentation, assert_class_issemiabstract, not_raises
from .common_fixtures_large import *

class DummyBinarySelectionProtocol(BinarySelectionProtocol):
    def problem(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
        pass

@pytest.fixture
def selsoln(
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
    out = DummyBinarySelectionProtocol(
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

################### Test class abstract/concrete properties ####################
def test_BinarySelectionProtocol_is_semiabstract():
    assert_class_issemiabstract(BinarySelectionProtocol)

############################## Test class docstring ############################
def test_BinarySelectionProtocol_docstring():
    assert_class_documentation(BinarySelectionProtocol)

############################## Test class methods ##############################
def test_BinarySelectionProtocol_problem_is_abstract():
    assert_method_isabstract(BinarySelectionProtocol, "problem")

def test_BinarySelectionProtocol_sosolve_is_concrete():
    assert_method_isconcrete(BinarySelectionProtocol, "sosolve")

def test_BinarySelectionProtocol_mosolve_is_concrete():
    assert_method_isconcrete(BinarySelectionProtocol, "mosolve")

def test_BinarySelectionProtocol_select_is_concrete():
    assert_method_isconcrete(BinarySelectionProtocol, "select")

############################# Test class utilities #############################
def test_check_is_BinarySelectionProtocol(selsoln):
    with not_raises(Exception):
        check_is_BinarySelectionProtocol(selsoln, "selsoln")
    with pytest.raises(TypeError):
        check_is_BinarySelectionProtocol(object(), "selsoln")
    with pytest.raises(TypeError):
        check_is_BinarySelectionProtocol(None, "selsoln")
