import numpy
import pytest

from pybrops.test.assert_python import assert_function_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_method_isabstract

from pybrops.breed.prot.sel.prob.BinarySelectionProblem import BinarySelectionProblem, check_is_BinarySelectionProblem

################################################################################
################################ Test fixtures #################################
################################################################################
class BinarySelectionProblemTestClass(BinarySelectionProblem):
    def __init__(
            self,
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, obj_trans, obj_trans_kwargs,
            nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs,
            neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs,
            **kwargs
        ):
        """NA"""
        super(BinarySelectionProblemTestClass, self).__init__(
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
        super(BinarySelectionProblemTestClass, self).latentfn(x, *args, **kwargs)

@pytest.fixture
def ndecn():
    yield 4

@pytest.fixture
def decn_space_lower():
    yield numpy.array([0,0,0,0], dtype=int)

@pytest.fixture
def decn_space_upper():
    yield numpy.array([1,1,1,1], dtype=int)

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
    yield BinarySelectionProblemTestClass(
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
def test_BinarySelectionProblem_docstring():
    assert_class_documentation(BinarySelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(BinarySelectionProblem, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_latentfn_is_abstract(prob):
    assert_method_isabstract(BinarySelectionProblem, "latentfn")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_BinarySelectionProblem_is_concrete():
    assert_function_isconcrete(check_is_BinarySelectionProblem)

def test_check_is_BinarySelectionProblem(prob):
    with not_raises(TypeError):
        check_is_BinarySelectionProblem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_BinarySelectionProblem(None, "prob")
