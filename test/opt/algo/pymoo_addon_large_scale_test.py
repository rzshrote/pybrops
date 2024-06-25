from typing import Iterable
import numpy
import pytest
from pybrops.opt.algo.pymoo_addon_large_scale import SubsetPoolReducedExchangeMutation, TiledSubsetPoolIterator
from pybrops.opt.algo.pymoo_addon_large_scale import SubsetPoolRandomSampling
from pybrops.test.assert_python import assert_class_isconcrete
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_module_documentation
from pybrops.test.assert_python import assert_module_public_api
from pybrops.test.assert_python import assert_property_isconcrete

################################################################################
################################# Mock Classes #################################
################################################################################

class DummyProblem:
    def __init__(self, n_var):
        self.n_var = n_var

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
################# TiledSubsetPoolIterator ##################
############################################################

@pytest.fixture
def nelem():
    yield 1000

@pytest.fixture
def setspace(nelem):
    out = numpy.arange(nelem)
    yield out

@pytest.fixture
def tspiter(setspace):
    out = TiledSubsetPoolIterator(setspace)
    yield out

############################################################
################ TiledSubsetRandomSampling #################
############################################################

@pytest.fixture
def tssample_n_var():
    yield 10

@pytest.fixture
def tssample_n_samples():
    yield 100

@pytest.fixture
def tssample_problem(tssample_n_var):
    out = DummyProblem(tssample_n_var)
    yield out

@pytest.fixture
def tssample(tspiter):
    out = SubsetPoolRandomSampling(tspiter)
    yield out

############################################################
############### TiledReducedExchangeMutation ###############
############################################################

@pytest.fixture
def tremutation(tspiter):
    out = SubsetPoolReducedExchangeMutation(tspiter)
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_pymoo_addon_large_scale_module_documentation():
    import pybrops.opt.algo.pymoo_addon_large_scale
    assert_module_documentation(pybrops.opt.algo.pymoo_addon_large_scale)

def test_pymoo_addon_large_scale_module_public_api():
    import pybrops.opt.algo.pymoo_addon_large_scale
    assert_module_public_api(pybrops.opt.algo.pymoo_addon_large_scale)

################################################################################
############################ Test class attributes #############################
################################################################################

def test_TiledChoicePoolIterator_is_concrete():
    assert_class_isconcrete(TiledSubsetPoolIterator)

def test_TiledSubsetRandomSampling_is_concrete():
    assert_class_isconcrete(SubsetPoolRandomSampling)

################################################################################
######################## Test abstract special methods #########################
################################################################################

################################################################################
######################## Test concrete special methods #########################
################################################################################

######
###### TiledSubsetPoolIterator
######

### __init__
def test_TiledChoicePoolIterator___init___is_concrete():
    assert_method_isconcrete(TiledSubsetPoolIterator, "__init__")

def test_TiledChoicePoolIterator___init__(tspiter):
    assert isinstance(tspiter, TiledSubsetPoolIterator)
    assert isinstance(tspiter._pool, numpy.ndarray)
    assert len(tspiter._pool) == len(tspiter.setspace)
    assert tspiter._poolix == 0

### __iter__
def test_TiledChoicePoolIterator___iter___is_concrete():
    assert_method_isconcrete(TiledSubsetPoolIterator, "__iter__")

def test_TiledChoicePoolIterator___iter__(tspiter):
    assert isinstance(iter(tspiter), Iterable)

### __next__
def test_TiledChoicePoolIterator___next___is_concrete():
    assert_method_isconcrete(TiledSubsetPoolIterator, "__next__")

def test_TiledChoicePoolIterator___next__(tspiter):
    assert numpy.issubdtype(next(tspiter), tspiter.setspace.dtype)

def test_TiledChoicePoolIterator___next___is_tiled(tspiter, nelem):
    observed = numpy.array([next(tspiter) for _ in range(nelem)])
    assert numpy.all(numpy.in1d(tspiter.setspace, observed))

######
###### TiledSubsetRandomSampling
######

### __init__
def test_TiledSubsetRandomSampling___init___is_concrete():
    assert_method_isconcrete(SubsetPoolRandomSampling, "__init__")

def test_TiledSubsetRandomSampling___init__(tssample):
    assert isinstance(tssample, SubsetPoolRandomSampling)

######
###### TiledReducedExchangeMutation
######

### __init__
def test_TiledReducedExchangeMutation___init___is_concrete():
    assert_method_isconcrete(SubsetPoolReducedExchangeMutation, "__init__")

def test_TiledReducedExchangeMutation___init__(tremutation):
    assert isinstance(tremutation, SubsetPoolReducedExchangeMutation)

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
########################### Test concrete properties ###########################
################################################################################

######
###### TiledSubsetPoolIterator
######

### setspace
def test_TiledChoicePoolIterator_setspace_is_concrete():
    assert_property_isconcrete(TiledSubsetPoolIterator, "setspace")

def test_TiledChoicePoolIterator_setspace(tspiter, setspace):
    assert numpy.all(tspiter.setspace == setspace)

######
###### TiledSubsetRandomSampling
######

### tspiter
def test_TiledSubsetRandomSampling_tspiter_is_concrete():
    assert_property_isconcrete(SubsetPoolRandomSampling, "tspiter")

######
###### TiledReducedExchangeMutation
######

### tspiter
def test_TiledReducedExchangeMutation_tspiter_is_concrete():
    assert_property_isconcrete(SubsetPoolReducedExchangeMutation, "tspiter")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################

######
###### TiledSubsetPoolIterator
######

### reset_pool
def test_TiledSubsetPoolIterator_reset_pool_is_concrete():
    assert_method_isconcrete(TiledSubsetPoolIterator, "reset_pool")

def test_TiledSubsetPoolIterator_reset_pool(tspiter):
    next(tspiter)
    old_pool = tspiter._pool.copy()
    old_poolix = tspiter._poolix
    tspiter.reset_pool()
    assert numpy.any(tspiter._pool != old_pool)
    assert tspiter._poolix != old_poolix
    assert tspiter._poolix == 0

### getnext
def test_TiledSubsetPoolIterator_getnext_is_concrete():
    assert_method_isconcrete(TiledSubsetPoolIterator, "getnext")

def test_TiledSubsetPoolIterator_getnext(tspiter):
    observed = tspiter.getnext(10)
    assert isinstance(observed, numpy.ndarray)
    assert len(observed) == 10

######
###### TiledSubsetRandomSampling
######

### _do
def test_TiledSubsetRandomSampling__do_is_concrete():
    assert_method_isconcrete(SubsetPoolRandomSampling, "_do")

def test_TiledSubsetRandomSampling__do(tssample, tssample_problem, tssample_n_samples, tssample_n_var):
    out = tssample._do(tssample_problem, tssample_n_samples)
    assert out.shape == (tssample_n_samples, tssample_n_var)

######
###### TiledReducedExchangeMutation
######

### _do
def test_TiledReducedExchangeMutation__do_is_concrete():
    assert_method_isconcrete(SubsetPoolReducedExchangeMutation, "_do")

def test_TiledReducedExchangeMutation__do(
        tssample, 
        tssample_problem, 
        tssample_n_samples, 
        tssample_n_var,
        tremutation
    ):
    X = tssample._do(tssample_problem, tssample_n_samples)
    Xm = tremutation._do(tssample_problem, X)
    assert Xm.shape == X.shape
    assert not numpy.all(X == Xm)

################################################################################
########################## Test abstract classmethods ##########################
################################################################################

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

################################################################################
######################### Test abstract staticmethods ##########################
################################################################################

################################################################################
######################### Test concrete staticmethods ##########################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
