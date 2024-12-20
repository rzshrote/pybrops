import numpy
import pytest

from pybrops.test.assert_python import assert_function_isconcrete, assert_class_documentation, not_raises
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_method_isabstract

from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem, check_is_SubsetSelectionProblem

################################################################################
################################ Test fixtures #################################
################################################################################
class SubsetSelectionProblemTestClass(SubsetSelectionProblem):
    def __init__(
            self,
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, obj_trans, obj_trans_kwargs,
            nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs,
            neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs,
            **kwargs
        ):
        """NA"""
        super(SubsetSelectionProblemTestClass, self).__init__(
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
        super(SubsetSelectionProblemTestClass, self).latentfn(x, *args, **kwargs)

@pytest.fixture
def ndecn():
    yield 4

@pytest.fixture
def decn_space_lower():
    yield numpy.array([1,2,3,4], dtype=float)

@pytest.fixture
def decn_space_upper():
    yield numpy.array([5,6,7,8], dtype=float)

@pytest.fixture
def decn_space(decn_space_lower, decn_space_upper):
    yield numpy.concatenate([decn_space_lower, decn_space_upper])

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
    yield SubsetSelectionProblemTestClass(
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
def test_SubsetSelectionProblem_docstring():
    assert_class_documentation(SubsetSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################
def test___init___is_concrete():
    assert_method_isconcrete(SubsetSelectionProblem, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_latentfn_is_abstract(prob):
    assert_method_isabstract(SubsetSelectionProblem, "latentfn")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SubsetSelectionProblem_is_concrete():
    assert_function_isconcrete(check_is_SubsetSelectionProblem)

def test_check_is_SubsetSelectionProblem(prob):
    with not_raises(TypeError):
        check_is_SubsetSelectionProblem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_SubsetSelectionProblem(None, "prob")
