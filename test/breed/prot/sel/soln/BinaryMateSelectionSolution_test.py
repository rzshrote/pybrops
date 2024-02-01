import pytest
from pybrops.breed.prot.sel.soln.BinaryMateSelectionSolution import BinaryMateSelectionSolution
from pybrops.breed.prot.sel.soln.BinaryMateSelectionSolution import BinaryMateSelectionSolution
from pybrops.breed.prot.sel.soln.BinaryMateSelectionSolution import check_is_BinaryMateSelectionSolution
from pybrops.test.assert_python import assert_docstring, assert_concrete_class, not_raises

from .common_fixtures import *

@pytest.fixture
def selsoln(
        common_ndecn_mate_binary,
        common_decn_space_mate_binary,
        common_decn_space_lower_mate_binary,
        common_decn_space_upper_mate_binary,
        common_decn_space_xmap,
        common_nobj,
        common_obj_wt,
        common_nineqcv,
        common_ineqcv_wt,
        common_neqcv,
        common_eqcv_wt,
        common_nsoln,
        common_soln_decn_mate_binary,
        common_soln_obj,
        common_soln_ineqcv,
        common_soln_eqcv
    ):
    out = BinaryMateSelectionSolution(
        ndecn = common_ndecn_mate_binary,
        decn_space = common_decn_space_mate_binary,
        decn_space_lower = common_decn_space_lower_mate_binary,
        decn_space_upper = common_decn_space_upper_mate_binary,
        decn_space_xmap = common_decn_space_xmap,
        nobj = common_nobj,
        obj_wt = common_obj_wt,
        nineqcv = common_nineqcv,
        ineqcv_wt = common_ineqcv_wt,
        neqcv = common_neqcv,
        eqcv_wt = common_eqcv_wt,
        nsoln = common_nsoln,
        soln_decn = common_soln_decn_mate_binary,
        soln_obj = common_soln_obj,
        soln_ineqcv = common_soln_ineqcv,
        soln_eqcv = common_soln_eqcv
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_BinaryMateSelectionSolution_is_concrete():
    assert_concrete_class(BinaryMateSelectionSolution)

############################## Test class docstring ############################
def test_BinaryMateSelectionSolution_docstring():
    assert_docstring(BinaryMateSelectionSolution)

############################# Test class utilities #############################
def test_check_is_BinaryMateSelectionSolution(selsoln):
    with not_raises(Exception):
        check_is_BinaryMateSelectionSolution(selsoln, "selsoln")
