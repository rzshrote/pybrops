import pytest
from pybrops.breed.prot.sel.soln.IntegerMateSelectionSolution import IntegerMateSelectionSolution
from pybrops.breed.prot.sel.soln.IntegerMateSelectionSolution import IntegerMateSelectionSolution
from pybrops.breed.prot.sel.soln.IntegerMateSelectionSolution import check_is_IntegerMateSelectionSolution
from pybrops.test.assert_python import assert_docstring, assert_concrete_class, not_raises

from pybrops.test.breed.prot.sel.soln.common_fixtures import *

@pytest.fixture
def selsoln(
        common_ndecn_mate_integer,
        common_decn_space_mate_integer,
        common_decn_space_lower_mate_integer,
        common_decn_space_upper_mate_integer,
        common_decn_space_xmap,
        common_nobj,
        common_obj_wt,
        common_nineqcv,
        common_ineqcv_wt,
        common_neqcv,
        common_eqcv_wt,
        common_nsoln,
        common_soln_decn_mate_integer,
        common_soln_obj,
        common_soln_ineqcv,
        common_soln_eqcv
    ):
    out = IntegerMateSelectionSolution(
        ndecn = common_ndecn_mate_integer,
        decn_space = common_decn_space_mate_integer,
        decn_space_lower = common_decn_space_lower_mate_integer,
        decn_space_upper = common_decn_space_upper_mate_integer,
        decn_space_xmap = common_decn_space_xmap,
        nobj = common_nobj,
        obj_wt = common_obj_wt,
        nineqcv = common_nineqcv,
        ineqcv_wt = common_ineqcv_wt,
        neqcv = common_neqcv,
        eqcv_wt = common_eqcv_wt,
        nsoln = common_nsoln,
        soln_decn = common_soln_decn_mate_integer,
        soln_obj = common_soln_obj,
        soln_ineqcv = common_soln_ineqcv,
        soln_eqcv = common_soln_eqcv
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_IntegerMateSelectionSolution_is_concrete():
    assert_concrete_class(IntegerMateSelectionSolution)

############################## Test class docstring ############################
def test_IntegerMateSelectionSolution_docstring():
    assert_docstring(IntegerMateSelectionSolution)

############################# Test class utilities #############################
def test_check_is_IntegerMateSelectionSolution(selsoln):
    with not_raises(Exception):
        check_is_IntegerMateSelectionSolution(selsoln, "selsoln")
