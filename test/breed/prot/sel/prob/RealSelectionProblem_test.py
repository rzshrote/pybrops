import numpy
import pytest

from pybrops.test.assert_python import assert_function_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_method_isabstract

from pybrops.breed.prot.sel.prob.RealSelectionProblem import RealSelectionProblem, check_is_RealSelectionProblem

################################################################################
################################ Test fixtures #################################
################################################################################
class RealSelectionProblemTestClass(RealSelectionProblem):
    def __init__(
            self,
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, obj_trans, obj_trans_kwargs,
            nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs,
            neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs,
            **kwargs
        ):
        """NA"""
        super(RealSelectionProblemTestClass, self).__init__(
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, obj_trans, obj_trans_kwargs,
            nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs,
            neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs,
            **kwargs
        )
    @property
    def nlatent(self) -> int:
        """nlatent."""
        return 0
    def latentfn(self, x, *args, **kwargs):
        """NA"""
        super(RealSelectionProblemTestClass, self).latentfn(x, *args, **kwargs)

@pytest.fixture
def ndecn():
    yield 4

@pytest.fixture
def decn_space_lower():
    yield numpy.array([0,1,2,3], dtype=float)

@pytest.fixture
def decn_space_upper():
    yield numpy.array([4,5,6,7], dtype=float)

@pytest.fixture
def decn_space(decn_space_lower, decn_space_upper):
    yield numpy.stack([decn_space_lower, decn_space_upper])

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
def prob(
    ndecn,
    decn_space,
    decn_space_lower,
    decn_space_upper,
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
    eqcv_trans_kwargs
):
    yield RealSelectionProblemTestClass(
        ndecn = ndecn,
        decn_space = decn_space,
        decn_space_lower = decn_space_lower,
        decn_space_upper = decn_space_upper,
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
        eqcv_trans_kwargs = eqcv_trans_kwargs
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_RealSelectionProblem_docstring():
    assert_class_documentation(RealSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################
def test___init___is_concrete():
    assert_method_isconcrete(RealSelectionProblem, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_latentfn_is_abstract(prob):
    assert_method_isabstract(RealSelectionProblem, "latentfn")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_RealSelectionProblem_is_concrete():
    assert_function_isconcrete(check_is_RealSelectionProblem)

def test_check_is_RealSelectionProblem(prob):
    with not_raises(TypeError):
        check_is_RealSelectionProblem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_RealSelectionProblem(None, "prob")
