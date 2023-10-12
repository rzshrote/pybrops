import numpy
import pytest
from pymoo.core.problem import ElementwiseEvaluationFunction, LoopedElementwiseEvaluation
from pybrops.core.random.prng import global_prng
from pybrops.opt.prob.BinaryProblem import BinaryProblem
from pybrops.opt.prob.IntegerProblem import IntegerProblem
from pybrops.opt.prob.RealProblem import RealProblem
from pybrops.opt.prob.SubsetProblem import SubsetProblem

@pytest.fixture
def common_ndecn():
    yield 20

@pytest.fixture
def common_superset_size(common_ndecn):
    yield common_ndecn * 3

@pytest.fixture
def common_decn_space_lower_binary(common_ndecn):
    yield numpy.repeat(0, common_ndecn)

@pytest.fixture
def common_decn_space_upper_binary(common_ndecn):
    yield numpy.repeat(1, common_ndecn)

@pytest.fixture
def common_decn_space_binary(common_decn_space_lower_binary, common_decn_space_upper_binary):
    yield numpy.stack([common_decn_space_lower_binary, common_decn_space_upper_binary])

@pytest.fixture
def common_decn_space_lower_integer(common_ndecn):
    yield numpy.repeat(0, common_ndecn)

@pytest.fixture
def common_decn_space_upper_integer(common_ndecn):
    yield numpy.repeat(common_ndecn, common_ndecn)

@pytest.fixture
def common_decn_space_integer(common_decn_space_lower_integer, common_decn_space_upper_integer):
    yield numpy.stack([common_decn_space_lower_integer, common_decn_space_upper_integer])

@pytest.fixture
def common_decn_space_lower_real(common_ndecn):
    yield numpy.repeat(0.0, common_ndecn)

@pytest.fixture
def common_decn_space_upper_real(common_ndecn):
    yield numpy.repeat(1.0, common_ndecn)

@pytest.fixture
def common_decn_space_real(common_decn_space_lower_real, common_decn_space_upper_real):
    yield numpy.stack([common_decn_space_lower_real, common_decn_space_upper_real])

@pytest.fixture
def common_decn_space_lower_subset(common_ndecn):
    yield numpy.repeat(0, common_ndecn)

@pytest.fixture
def common_decn_space_upper_subset(common_ndecn, common_superset_size):
    yield numpy.repeat(common_superset_size-1, common_ndecn)

@pytest.fixture
def common_decn_space_subset(common_superset_size):
    yield numpy.arange(common_superset_size)

@pytest.fixture
def common_nobj_single():
    yield 1

@pytest.fixture
def common_obj_wt_single():
    yield numpy.array([1], dtype=float)

@pytest.fixture
def common_nineqcv_single():
    yield 0

@pytest.fixture
def common_ineqcv_wt_single():
    yield numpy.array([], dtype=float)

@pytest.fixture
def common_neqcv_single():
    yield 0

@pytest.fixture
def common_eqcv_wt_single():
    yield numpy.array([], dtype=float)

@pytest.fixture
def common_nobj_multi():
    yield 2

@pytest.fixture
def common_obj_wt_multi():
    yield numpy.array([1, 1], dtype=float)

@pytest.fixture
def common_nineqcv_multi():
    yield 0

@pytest.fixture
def common_ineqcv_wt_multi():
    yield numpy.array([], dtype=float)

@pytest.fixture
def common_neqcv_multi():
    yield 0

@pytest.fixture
def common_eqcv_wt_multi():
    yield numpy.array([], dtype=float)

### problems ###
class DummySingleObjectiveSummationBinaryProblem(BinaryProblem):
    def __init__(
            self, ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        ):
        """NA"""
        super(DummySingleObjectiveSummationBinaryProblem, self).__init__(
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        )
    ### method required by PyBrOpS interface ###
    def evalfn(self, x, *args, **kwargs):
        """NA"""
        obj = self.obj_wt * numpy.sum(x, axis=0, keepdims=True)
        ineqcv = self.ineqcv_wt * numpy.zeros(self.ineqcv_wt)
        eqcv = self.eqcv_wt * numpy.zeros(self.eqcv_wt)
        # return results
        return obj, ineqcv, eqcv
    ### method required by PyMOO interface ###
    def _evaluate(self, x, out, *args, **kwargs):
        """NA"""
        # if x is a vector, score and update output dictionary
        if x.ndim == 1:
            vals = self.evalfn(x, *args, **kwargs)
            out.update({key:val for key,val in zip(["F","G","H"],vals) if len(val) > 0})
        # if x is a matrix or other, score each row and update output dictionary
        else:
            vals = [self.evalfn(v *args, **kwargs) for v in x]  # evaluate each vector
            obj = numpy.stack([e[0] for e in vals])             # extract objective function evaluations
            ineqcv = numpy.stack([e[1] for e in vals])          # extract inequality constraint evaluations
            eqcv = numpy.stack([e[2] for e in vals])            # extract equality constraint evaluations
            out.update({key:val for key,val in zip(["F","G","H"],[obj,ineqcv,eqcv]) if val.shape[1] > 0})

class DummyMultiObjectiveSummationBinaryProblem(BinaryProblem):
    def __init__(
            self, ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        ):
        """NA"""
        super(DummyMultiObjectiveSummationBinaryProblem, self).__init__(
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        )
    ### method required by PyBrOpS interface ###
    def evalfn(self, x, *args, **kwargs):
        """NA"""
        obj = self.obj_wt * numpy.array([numpy.sum(x), numpy.sum(1-x)], dtype=float)
        ineqcv = self.ineqcv_wt * numpy.zeros(self.ineqcv_wt)
        eqcv = self.eqcv_wt * numpy.zeros(self.eqcv_wt)
        # return results
        return obj, ineqcv, eqcv
    ### method required by PyMOO interface ###
    def _evaluate(self, x, out, *args, **kwargs):
        """NA"""
        # if x is a vector, score and update output dictionary
        if x.ndim == 1:
            vals = self.evalfn(x, *args, **kwargs)
            out.update({key:val for key,val in zip(["F","G","H"],vals) if len(val) > 0})
        # if x is a matrix or other, score each row and update output dictionary
        else:
            vals = [self.evalfn(v *args, **kwargs) for v in x]  # evaluate each vector
            obj = numpy.stack([e[0] for e in vals])             # extract objective function evaluations
            ineqcv = numpy.stack([e[1] for e in vals])          # extract inequality constraint evaluations
            eqcv = numpy.stack([e[2] for e in vals])            # extract equality constraint evaluations
            out.update({key:val for key,val in zip(["F","G","H"],[obj,ineqcv,eqcv]) if val.shape[1] > 0})

@pytest.fixture
def common_sum_prob_binary_single(
        common_ndecn,
        common_decn_space_binary,
        common_decn_space_lower_binary,
        common_decn_space_upper_binary,
        common_nobj_single,
        common_obj_wt_single,
        common_nineqcv_single,
        common_ineqcv_wt_single,
        common_neqcv_single,
        common_eqcv_wt_single
    ):
    out = DummySingleObjectiveSummationBinaryProblem(
        ndecn = common_ndecn,
        decn_space = common_decn_space_binary,
        decn_space_lower = common_decn_space_lower_binary,
        decn_space_upper = common_decn_space_upper_binary,
        nobj = common_nobj_single,
        obj_wt = common_obj_wt_single,
        nineqcv = common_nineqcv_single,
        ineqcv_wt = common_ineqcv_wt_single,
        neqcv = common_neqcv_single,
        eqcv_wt = common_eqcv_wt_single,
        vtype = None,
        vars = None,
        elementwise = True,
        elementwise_func = ElementwiseEvaluationFunction,
        elementwise_runner = LoopedElementwiseEvaluation(),
        replace_nan_values_by = None,
        exclude_from_serialization = None,
        callback = None,
        strict = True,
    )
    yield out

@pytest.fixture
def common_sum_prob_binary_multi(
        common_ndecn,
        common_decn_space_binary,
        common_decn_space_lower_binary,
        common_decn_space_upper_binary,
        common_nobj_multi,
        common_obj_wt_multi,
        common_nineqcv_multi,
        common_ineqcv_wt_multi,
        common_neqcv_multi,
        common_eqcv_wt_multi
    ):
    out = DummyMultiObjectiveSummationBinaryProblem(
        ndecn = common_ndecn,
        decn_space = common_decn_space_binary,
        decn_space_lower = common_decn_space_lower_binary,
        decn_space_upper = common_decn_space_upper_binary,
        nobj = common_nobj_multi,
        obj_wt = common_obj_wt_multi,
        nineqcv = common_nineqcv_multi,
        ineqcv_wt = common_ineqcv_wt_multi,
        neqcv = common_neqcv_multi,
        eqcv_wt = common_eqcv_wt_multi,
        vtype = None,
        vars = None,
        elementwise = True,
        elementwise_func = ElementwiseEvaluationFunction,
        elementwise_runner = LoopedElementwiseEvaluation(),
        replace_nan_values_by = None,
        exclude_from_serialization = None,
        callback = None,
        strict = True,
    )
    yield out

class DummySingleObjectiveSummationIntegerProblem(IntegerProblem):
    def __init__(
            self, ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        ):
        """NA"""
        super(DummySingleObjectiveSummationIntegerProblem, self).__init__(
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        )
    ### method required by PyBrOpS interface ###
    def evalfn(self, x, *args, **kwargs):
        """NA"""
        obj = self.obj_wt * numpy.sum(x, axis=0, keepdims=True)
        ineqcv = self.ineqcv_wt * numpy.zeros(self.ineqcv_wt)
        eqcv = self.eqcv_wt * numpy.zeros(self.eqcv_wt)
        # return results
        return obj, ineqcv, eqcv
    ### method required by PyMOO interface ###
    def _evaluate(self, x, out, *args, **kwargs):
        """NA"""
        # if x is a vector, score and update output dictionary
        if x.ndim == 1:
            vals = self.evalfn(x, *args, **kwargs)
            out.update({key:val for key,val in zip(["F","G","H"],vals) if len(val) > 0})
        # if x is a matrix or other, score each row and update output dictionary
        else:
            vals = [self.evalfn(v *args, **kwargs) for v in x]  # evaluate each vector
            obj = numpy.stack([e[0] for e in vals])             # extract objective function evaluations
            ineqcv = numpy.stack([e[1] for e in vals])          # extract inequality constraint evaluations
            eqcv = numpy.stack([e[2] for e in vals])            # extract equality constraint evaluations
            out.update({key:val for key,val in zip(["F","G","H"],[obj,ineqcv,eqcv]) if val.shape[1] > 0})

class DummyMultiObjectiveSummationIntegerProblem(IntegerProblem):
    def __init__(
            self, ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        ):
        """NA"""
        super(DummyMultiObjectiveSummationIntegerProblem, self).__init__(
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        )
    ### method required by PyBrOpS interface ###
    def evalfn(self, x, *args, **kwargs):
        """NA"""
        obj = self.obj_wt * numpy.array([numpy.sum(x), numpy.sum(1-x)], dtype=float)
        ineqcv = self.ineqcv_wt * numpy.zeros(self.ineqcv_wt)
        eqcv = self.eqcv_wt * numpy.zeros(self.eqcv_wt)
        # return results
        return obj, ineqcv, eqcv
    ### method required by PyMOO interface ###
    def _evaluate(self, x, out, *args, **kwargs):
        """NA"""
        # if x is a vector, score and update output dictionary
        if x.ndim == 1:
            vals = self.evalfn(x, *args, **kwargs)
            out.update({key:val for key,val in zip(["F","G","H"],vals) if len(val) > 0})
        # if x is a matrix or other, score each row and update output dictionary
        else:
            vals = [self.evalfn(v *args, **kwargs) for v in x]  # evaluate each vector
            obj = numpy.stack([e[0] for e in vals])             # extract objective function evaluations
            ineqcv = numpy.stack([e[1] for e in vals])          # extract inequality constraint evaluations
            eqcv = numpy.stack([e[2] for e in vals])            # extract equality constraint evaluations
            out.update({key:val for key,val in zip(["F","G","H"],[obj,ineqcv,eqcv]) if val.shape[1] > 0})

@pytest.fixture
def common_sum_prob_integer_single(
        common_ndecn,
        common_decn_space_integer,
        common_decn_space_lower_integer,
        common_decn_space_upper_integer,
        common_nobj_single,
        common_obj_wt_single,
        common_nineqcv_single,
        common_ineqcv_wt_single,
        common_neqcv_single,
        common_eqcv_wt_single
    ):
    out = DummySingleObjectiveSummationIntegerProblem(
        ndecn = common_ndecn,
        decn_space = common_decn_space_integer,
        decn_space_lower = common_decn_space_lower_integer,
        decn_space_upper = common_decn_space_upper_integer,
        nobj = common_nobj_single,
        obj_wt = common_obj_wt_single,
        nineqcv = common_nineqcv_single,
        ineqcv_wt = common_ineqcv_wt_single,
        neqcv = common_neqcv_single,
        eqcv_wt = common_eqcv_wt_single,
        vtype = None,
        vars = None,
        elementwise = True,
        elementwise_func = ElementwiseEvaluationFunction,
        elementwise_runner = LoopedElementwiseEvaluation(),
        replace_nan_values_by = None,
        exclude_from_serialization = None,
        callback = None,
        strict = True,
    )
    yield out

@pytest.fixture
def common_sum_prob_integer_multi(
        common_ndecn,
        common_decn_space_integer,
        common_decn_space_lower_integer,
        common_decn_space_upper_integer,
        common_nobj_multi,
        common_obj_wt_multi,
        common_nineqcv_multi,
        common_ineqcv_wt_multi,
        common_neqcv_multi,
        common_eqcv_wt_multi
    ):
    out = DummyMultiObjectiveSummationIntegerProblem(
        ndecn = common_ndecn,
        decn_space = common_decn_space_integer,
        decn_space_lower = common_decn_space_lower_integer,
        decn_space_upper = common_decn_space_upper_integer,
        nobj = common_nobj_multi,
        obj_wt = common_obj_wt_multi,
        nineqcv = common_nineqcv_multi,
        ineqcv_wt = common_ineqcv_wt_multi,
        neqcv = common_neqcv_multi,
        eqcv_wt = common_eqcv_wt_multi,
        vtype = None,
        vars = None,
        elementwise = True,
        elementwise_func = ElementwiseEvaluationFunction,
        elementwise_runner = LoopedElementwiseEvaluation(),
        replace_nan_values_by = None,
        exclude_from_serialization = None,
        callback = None,
        strict = True,
    )
    yield out

class DummySingleObjectiveSummationRealProblem(RealProblem):
    def __init__(
            self, ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        ):
        """NA"""
        super(DummySingleObjectiveSummationRealProblem, self).__init__(
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        )
    ### method required by PyBrOpS interface ###
    def evalfn(self, x, *args, **kwargs):
        """NA"""
        obj = self.obj_wt * numpy.sum(x, axis=0, keepdims=True)
        ineqcv = self.ineqcv_wt * numpy.zeros(self.ineqcv_wt)
        eqcv = self.eqcv_wt * numpy.zeros(self.eqcv_wt)
        # return results
        return obj, ineqcv, eqcv
    ### method required by PyMOO interface ###
    def _evaluate(self, x, out, *args, **kwargs):
        """NA"""
        # if x is a vector, score and update output dictionary
        if x.ndim == 1:
            vals = self.evalfn(x, *args, **kwargs)
            out.update({key:val for key,val in zip(["F","G","H"],vals) if len(val) > 0})
        # if x is a matrix or other, score each row and update output dictionary
        else:
            vals = [self.evalfn(v *args, **kwargs) for v in x]  # evaluate each vector
            obj = numpy.stack([e[0] for e in vals])             # extract objective function evaluations
            ineqcv = numpy.stack([e[1] for e in vals])          # extract inequality constraint evaluations
            eqcv = numpy.stack([e[2] for e in vals])            # extract equality constraint evaluations
            out.update({key:val for key,val in zip(["F","G","H"],[obj,ineqcv,eqcv]) if val.shape[1] > 0})

class DummyMultiObjectiveSummationRealProblem(RealProblem):
    def __init__(
            self, ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        ):
        """NA"""
        super(DummyMultiObjectiveSummationRealProblem, self).__init__(
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        )
    ### method required by PyBrOpS interface ###
    def evalfn(self, x, *args, **kwargs):
        """NA"""
        obj = self.obj_wt * numpy.array([numpy.sum(x), numpy.sum(1-x)], dtype=float)
        ineqcv = self.ineqcv_wt * numpy.zeros(self.ineqcv_wt)
        eqcv = self.eqcv_wt * numpy.zeros(self.eqcv_wt)
        # return results
        return obj, ineqcv, eqcv
    ### method required by PyMOO interface ###
    def _evaluate(self, x, out, *args, **kwargs):
        """NA"""
        # if x is a vector, score and update output dictionary
        if x.ndim == 1:
            vals = self.evalfn(x, *args, **kwargs)
            out.update({key:val for key,val in zip(["F","G","H"],vals) if len(val) > 0})
        # if x is a matrix or other, score each row and update output dictionary
        else:
            vals = [self.evalfn(v *args, **kwargs) for v in x]  # evaluate each vector
            obj = numpy.stack([e[0] for e in vals])             # extract objective function evaluations
            ineqcv = numpy.stack([e[1] for e in vals])          # extract inequality constraint evaluations
            eqcv = numpy.stack([e[2] for e in vals])            # extract equality constraint evaluations
            out.update({key:val for key,val in zip(["F","G","H"],[obj,ineqcv,eqcv]) if val.shape[1] > 0})

@pytest.fixture
def common_sum_prob_real_single(
        common_ndecn,
        common_decn_space_real,
        common_decn_space_lower_real,
        common_decn_space_upper_real,
        common_nobj_single,
        common_obj_wt_single,
        common_nineqcv_single,
        common_ineqcv_wt_single,
        common_neqcv_single,
        common_eqcv_wt_single
    ):
    out = DummySingleObjectiveSummationRealProblem(
        ndecn = common_ndecn,
        decn_space = common_decn_space_real,
        decn_space_lower = common_decn_space_lower_real,
        decn_space_upper = common_decn_space_upper_real,
        nobj = common_nobj_single,
        obj_wt = common_obj_wt_single,
        nineqcv = common_nineqcv_single,
        ineqcv_wt = common_ineqcv_wt_single,
        neqcv = common_neqcv_single,
        eqcv_wt = common_eqcv_wt_single,
        vtype = None,
        vars = None,
        elementwise = True,
        elementwise_func = ElementwiseEvaluationFunction,
        elementwise_runner = LoopedElementwiseEvaluation(),
        replace_nan_values_by = None,
        exclude_from_serialization = None,
        callback = None,
        strict = True,
    )
    yield out

@pytest.fixture
def common_sum_prob_real_multi(
        common_ndecn,
        common_decn_space_real,
        common_decn_space_lower_real,
        common_decn_space_upper_real,
        common_nobj_multi,
        common_obj_wt_multi,
        common_nineqcv_multi,
        common_ineqcv_wt_multi,
        common_neqcv_multi,
        common_eqcv_wt_multi
    ):
    out = DummyMultiObjectiveSummationRealProblem(
        ndecn = common_ndecn,
        decn_space = common_decn_space_real,
        decn_space_lower = common_decn_space_lower_real,
        decn_space_upper = common_decn_space_upper_real,
        nobj = common_nobj_multi,
        obj_wt = common_obj_wt_multi,
        nineqcv = common_nineqcv_multi,
        ineqcv_wt = common_ineqcv_wt_multi,
        neqcv = common_neqcv_multi,
        eqcv_wt = common_eqcv_wt_multi,
        vtype = None,
        vars = None,
        elementwise = True,
        elementwise_func = ElementwiseEvaluationFunction,
        elementwise_runner = LoopedElementwiseEvaluation(),
        replace_nan_values_by = None,
        exclude_from_serialization = None,
        callback = None,
        strict = True,
    )
    yield out

class DummySingleObjectiveSummationSubsetProblem(SubsetProblem):
    def __init__(
            self, ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        ):
        """NA"""
        super(DummySingleObjectiveSummationSubsetProblem, self).__init__(
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        )
    ### method required by PyBrOpS interface ###
    def evalfn(self, x, *args, **kwargs):
        """NA"""
        obj = self.obj_wt * numpy.sum(x, axis=0, keepdims=True)
        ineqcv = self.ineqcv_wt * numpy.zeros(self.ineqcv_wt)
        eqcv = self.eqcv_wt * numpy.zeros(self.eqcv_wt)
        # return results
        return obj, ineqcv, eqcv
    ### method required by PyMOO interface ###
    def _evaluate(self, x, out, *args, **kwargs):
        """NA"""
        # if x is a vector, score and update output dictionary
        if x.ndim == 1:
            vals = self.evalfn(x, *args, **kwargs)
            out.update({key:val for key,val in zip(["F","G","H"],vals) if len(val) > 0})
        # if x is a matrix or other, score each row and update output dictionary
        else:
            vals = [self.evalfn(v *args, **kwargs) for v in x]  # evaluate each vector
            obj = numpy.stack([e[0] for e in vals])             # extract objective function evaluations
            ineqcv = numpy.stack([e[1] for e in vals])          # extract inequality constraint evaluations
            eqcv = numpy.stack([e[2] for e in vals])            # extract equality constraint evaluations
            out.update({key:val for key,val in zip(["F","G","H"],[obj,ineqcv,eqcv]) if val.shape[1] > 0})

class DummyMultiObjectiveSummationSubsetProblem(SubsetProblem):
    def __init__(
            self, ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        ):
        """NA"""
        super(DummyMultiObjectiveSummationSubsetProblem, self).__init__(
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        )
    ### method required by PyBrOpS interface ###
    def evalfn(self, x, *args, **kwargs):
        """NA"""
        obj = self.obj_wt * numpy.array([numpy.sum(x), numpy.sum(1-x)], dtype=float)
        ineqcv = self.ineqcv_wt * numpy.zeros(self.ineqcv_wt)
        eqcv = self.eqcv_wt * numpy.zeros(self.eqcv_wt)
        # return results
        return obj, ineqcv, eqcv
    ### method required by PyMOO interface ###
    def _evaluate(self, x, out, *args, **kwargs):
        """NA"""
        # if x is a vector, score and update output dictionary
        if x.ndim == 1:
            vals = self.evalfn(x, *args, **kwargs)
            out.update({key:val for key,val in zip(["F","G","H"],vals) if len(val) > 0})
        # if x is a matrix or other, score each row and update output dictionary
        else:
            vals = [self.evalfn(v *args, **kwargs) for v in x]  # evaluate each vector
            obj = numpy.stack([e[0] for e in vals])             # extract objective function evaluations
            ineqcv = numpy.stack([e[1] for e in vals])          # extract inequality constraint evaluations
            eqcv = numpy.stack([e[2] for e in vals])            # extract equality constraint evaluations
            out.update({key:val for key,val in zip(["F","G","H"],[obj,ineqcv,eqcv]) if val.shape[1] > 0})

@pytest.fixture
def common_sum_prob_subset_single(
        common_ndecn,
        common_decn_space_subset,
        common_decn_space_lower_subset,
        common_decn_space_upper_subset,
        common_nobj_single,
        common_obj_wt_single,
        common_nineqcv_single,
        common_ineqcv_wt_single,
        common_neqcv_single,
        common_eqcv_wt_single
    ):
    out = DummySingleObjectiveSummationSubsetProblem(
        ndecn = common_ndecn,
        decn_space = common_decn_space_subset,
        decn_space_lower = common_decn_space_lower_subset,
        decn_space_upper = common_decn_space_upper_subset,
        nobj = common_nobj_single,
        obj_wt = common_obj_wt_single,
        nineqcv = common_nineqcv_single,
        ineqcv_wt = common_ineqcv_wt_single,
        neqcv = common_neqcv_single,
        eqcv_wt = common_eqcv_wt_single,
        vtype = None,
        vars = None,
        elementwise = True,
        elementwise_func = ElementwiseEvaluationFunction,
        elementwise_runner = LoopedElementwiseEvaluation(),
        replace_nan_values_by = None,
        exclude_from_serialization = None,
        callback = None,
        strict = True,
    )
    yield out

@pytest.fixture
def common_sum_prob_subset_multi(
        common_ndecn,
        common_decn_space_subset,
        common_decn_space_lower_subset,
        common_decn_space_upper_subset,
        common_nobj_multi,
        common_obj_wt_multi,
        common_nineqcv_multi,
        common_ineqcv_wt_multi,
        common_neqcv_multi,
        common_eqcv_wt_multi
    ):
    out = DummyMultiObjectiveSummationSubsetProblem(
        ndecn = common_ndecn,
        decn_space = common_decn_space_subset,
        decn_space_lower = common_decn_space_lower_subset,
        decn_space_upper = common_decn_space_upper_subset,
        nobj = common_nobj_multi,
        obj_wt = common_obj_wt_multi,
        nineqcv = common_nineqcv_multi,
        ineqcv_wt = common_ineqcv_wt_multi,
        neqcv = common_neqcv_multi,
        eqcv_wt = common_eqcv_wt_multi,
        vtype = None,
        vars = None,
        elementwise = True,
        elementwise_func = ElementwiseEvaluationFunction,
        elementwise_runner = LoopedElementwiseEvaluation(),
        replace_nan_values_by = None,
        exclude_from_serialization = None,
        callback = None,
        strict = True,
    )
    yield out

### algorithm parameters ###
@pytest.fixture
def common_ngen():
    yield 250

@pytest.fixture
def common_pop_size():
    yield 100

@pytest.fixture
def common_rng():
    yield global_prng
