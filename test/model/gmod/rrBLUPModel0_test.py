import copy
from numbers import Real
import os
import numpy
import pandas
import pytest
from matplotlib import pyplot
from pybrops.model.gmod.rrBLUPModel0 import check_is_rrBLUPModel0, rrBLUP_ML0, rrBLUPModel0
from pybrops.model.gmod.rrBLUPModel0 import rrBLUP_ML0_calc_G
from pybrops.model.gmod.rrBLUPModel0 import rrBLUP_ML0_calc_d_V
from pybrops.model.gmod.rrBLUPModel0 import rrBLUP_ML0_calc_etasq
from pybrops.model.gmod.rrBLUPModel0 import rrBLUP_ML0_center_y
from pybrops.model.gmod.rrBLUPModel0 import rrBLUP_ML0_neg2LogLik_fast
from pybrops.model.gmod.rrBLUPModel0 import rrBLUP_ML0_nonzero_d_V
from pybrops.model.gmod.rrBLUPModel0 import rrBLUP_ML0_calc_ridge
from pybrops.model.gmod.rrBLUPModel0 import rrBLUP_ML0_calc_ZtZplI
from pybrops.model.gmod.rrBLUPModel0 import rrBLUP_ML0_calc_Zty
from pybrops.model.gmod.rrBLUPModel0 import gauss_seidel
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.test.assert_python import assert_concrete_function, assert_concrete_method, assert_docstring, not_raises

################################ Test fixtures #################################

###################### Wheat Fixtures ######################

@pytest.fixture
def wheat_yvec():
    traits_df = pandas.read_csv("data/traits.csv")
    out = traits_df.iloc[:,1].to_numpy(dtype = float)
    yield out

@pytest.fixture
def wheat_Ymat(wheat_yvec):
    out = wheat_yvec[:,None]
    yield out

@pytest.fixture
def wheat_Zmat():
    markers_df = pandas.read_csv("data/markers.csv")
    out = markers_df.to_numpy(dtype = float)
    yield out

@pytest.fixture
def wheat_nobs(wheat_yvec):
    yield len(wheat_yvec)

@pytest.fixture
def wheat_nmkr(wheat_Zmat):
    yield wheat_Zmat.shape[1]

@pytest.fixture
def wheat_varE():
    yield 0.5725375

@pytest.fixture
def wheat_varU():
    yield 0.4690475

################# General shape parameters #################
@pytest.fixture
def ntrait():
    yield 3

@pytest.fixture
def ntaxa():
    yield 100

@pytest.fixture
def nmarker():
    yield 100

@pytest.fixture
def nphase():
    yield 2

@pytest.fixture
def nchrom():
    yield 10

@pytest.fixture
def ngroup():
    yield 5

###################### Genomic model #######################
@pytest.fixture
def algmod_beta(ntrait):
    out = numpy.random.uniform(0, 10, (1,ntrait))
    yield out

@pytest.fixture
def algmod_u_misc():
    yield None

@pytest.fixture
def algmod_u_a(nmarker, ntrait):
    out = numpy.random.normal(size = (nmarker, ntrait))
    yield out

@pytest.fixture
def algmod_trait(ntrait):
    out = numpy.array(["Trait"+str(i+1) for i in range(ntrait)], dtype=object)
    yield out

@pytest.fixture
def algmod_model_name():
    yield "test_algmod"

@pytest.fixture
def algmod_params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def algmod(
        algmod_beta, 
        algmod_u_misc, 
        algmod_u_a, 
        algmod_trait, 
        algmod_model_name, 
        algmod_params
    ):
    yield rrBLUPModel0(
        beta       = algmod_beta,
        u_misc     = algmod_u_misc,
        u_a        = algmod_u_a,
        trait      = algmod_trait,
        model_name = algmod_model_name,
        hyperparams     = algmod_params
    )

######################## Genotypes #########################
@pytest.fixture
def pgmat_mat(nphase, ntaxa, nmarker):
    out = numpy.random.randint(0,2,(nphase,ntaxa,nmarker))
    out = out.astype("int8")
    yield out

@pytest.fixture
def pgmat_chrgrp(nchrom, nmarker):
    out = numpy.random.randint(1, nchrom+1, nmarker)
    out.sort()
    yield out

@pytest.fixture
def pgmat_phypos(nmarker):
    out = numpy.random.randint(1, 2**32-1, nmarker)
    out.sort()
    yield out

@pytest.fixture
def pgmat_taxa(ntaxa):
    out = numpy.array(["Line"+str(i+1).zfill(3) for i in range(ntaxa)], dtype=object)
    yield out

@pytest.fixture
def pgmat_taxa_grp(ngroup, ntaxa):
    out = numpy.random.randint(1, ngroup+1, ntaxa)
    out.sort()
    yield out

@pytest.fixture
def pgmat(
        pgmat_mat, 
        pgmat_chrgrp, 
        pgmat_phypos, 
        pgmat_taxa, 
        pgmat_taxa_grp
    ):
    yield DensePhasedGenotypeMatrix(
        mat = pgmat_mat,
        vrnt_chrgrp = pgmat_chrgrp,
        vrnt_phypos = pgmat_phypos,
        taxa = pgmat_taxa,
        taxa_grp = pgmat_taxa_grp
    )

##################### Breeding values ######################
@pytest.fixture
def mat_intercept(pgmat, algmod_beta):
    n = pgmat.ntaxa
    q = algmod_beta.shape[0]
    a = numpy.ones((n,1), dtype = "float64")
    yield a

@pytest.fixture
def bvmat(algmod, pgmat):
    out = algmod.gebv(pgmat)
    sd = out.mat.std(0)
    sd *= 0.5
    out.mat += numpy.random.normal(0, sd, (out.ntaxa,len(sd)))
    yield out

#############################
### Helper Function Tests ###
#############################

def test_rrBLUP_ML0_calc_G(wheat_Zmat, wheat_nobs):
    G = rrBLUP_ML0_calc_G(wheat_Zmat)
    assert isinstance(G, numpy.ndarray)
    assert G.shape == (wheat_nobs,wheat_nobs)

def test_rrBLUP_ML0_center_y(wheat_yvec, wheat_nobs):
    y = rrBLUP_ML0_center_y(wheat_yvec)
    assert isinstance(y, numpy.ndarray)
    assert y.shape == (wheat_nobs,)
    assert abs(y.mean(0)) < 1e-10

def test_rrBLUP_ML0_calc_d_V(wheat_Zmat, wheat_nobs):
    G = rrBLUP_ML0_calc_G(wheat_Zmat)
    out = rrBLUP_ML0_calc_d_V(G)
    assert isinstance(out, tuple)
    assert len(out) == 2
    d = out[0]
    assert isinstance(d, numpy.ndarray)
    assert d.shape == (wheat_nobs,)
    V = out[1]
    assert isinstance(V, numpy.ndarray)
    assert V.shape == (wheat_nobs,wheat_nobs)

def test_rrBLUP_ML0_nonzero_d_V(wheat_Zmat, wheat_nobs):
    G = rrBLUP_ML0_calc_G(wheat_Zmat)
    d, V = rrBLUP_ML0_calc_d_V(G)
    out = rrBLUP_ML0_nonzero_d_V(d, V)
    assert isinstance(out, tuple)
    assert len(out) == 2
    d = out[0]
    assert isinstance(d, numpy.ndarray)
    assert d.shape == (wheat_nobs-1,)
    V = out[1]
    assert isinstance(V, numpy.ndarray)
    assert V.shape == (wheat_nobs,wheat_nobs-1)

def test_rrBLUP_ML0_calc_etasq(wheat_yvec, wheat_Zmat, wheat_nobs):
    G = rrBLUP_ML0_calc_G(wheat_Zmat)
    d, V = rrBLUP_ML0_calc_d_V(G)
    d, V = rrBLUP_ML0_nonzero_d_V(d, V)
    y = rrBLUP_ML0_center_y(wheat_yvec)
    etasq = rrBLUP_ML0_calc_etasq(y, V)
    assert isinstance(etasq, numpy.ndarray)
    assert etasq.shape == (wheat_nobs-1,)

def test_rrBLUP_ML0_neg2logLik_fast(wheat_yvec, wheat_Zmat, wheat_nobs, wheat_varE, wheat_varU):
    G = rrBLUP_ML0_calc_G(wheat_Zmat)
    d, V = rrBLUP_ML0_calc_d_V(G)
    d, V = rrBLUP_ML0_nonzero_d_V(d, V)
    y = rrBLUP_ML0_center_y(wheat_yvec)
    etasq = rrBLUP_ML0_calc_etasq(y, V)
    out = rrBLUP_ML0_neg2LogLik_fast(numpy.log([wheat_varE, wheat_varU]), etasq, d, wheat_nobs)
    assert isinstance(out, Real)
    # plot contour of function for visual examination
    pts = numpy.linspace(-1, 0, 30)
    gridpts = numpy.meshgrid(pts, pts)
    gridX = gridpts[0] # (g,g) containing log(wheat_varE) values
    gridY = gridpts[1] # (g,g) containing log(wheat_varU) values
    gridZ = numpy.empty(gridX.shape, dtype = float)
    for i in range(gridX.shape[0]):
        for j in range(gridX.shape[1]):
            gridZ[i,j] = rrBLUP_ML0_neg2LogLik_fast((gridX[i,j],gridY[i,j]), etasq, d, len(y))
    fig, ax = pyplot.subplots()
    CS = ax.contour(gridX, gridY, gridZ, levels = 10)
    ax.clabel(CS, inline=True, fontsize=10)
    ax.set_xlabel("log(wheat_varE)")
    ax.set_ylabel("log(wheat_varU)")
    ax.set_title('-2 * log-likelihood (minimizing)')
    pyplot.savefig("neg2LogLik.png")
    pyplot.close()
    
def test_rrBLUP_ML0_calc_ridge(wheat_varE, wheat_varU):
    ridge = rrBLUP_ML0_calc_ridge(wheat_varE, wheat_varU)
    assert isinstance(ridge, Real)
    assert ridge == (wheat_varE / wheat_varU)

def test_rrBLUP_ML0_calc_ZtZplI(wheat_Zmat, wheat_varE, wheat_varU, wheat_nmkr):
    ridge = rrBLUP_ML0_calc_ridge(wheat_varE, wheat_varU)
    A = rrBLUP_ML0_calc_ZtZplI(wheat_Zmat, ridge)
    assert isinstance(A, numpy.ndarray)
    assert A.shape == (wheat_nmkr,wheat_nmkr)

def test_rrBLUP_ML0_calc_Zty(wheat_Zmat, wheat_yvec, wheat_nmkr):
    y = rrBLUP_ML0_center_y(wheat_yvec)
    b = rrBLUP_ML0_calc_Zty(wheat_Zmat, y)
    assert isinstance(b, numpy.ndarray)
    assert b.shape == (wheat_nmkr,)

def test_gauss_seidel(wheat_yvec, wheat_Zmat, wheat_nmkr, wheat_varE, wheat_varU):
    ridge = rrBLUP_ML0_calc_ridge(wheat_varE, wheat_varU)
    A = rrBLUP_ML0_calc_ZtZplI(wheat_Zmat, ridge)
    y = rrBLUP_ML0_center_y(wheat_yvec)
    b = rrBLUP_ML0_calc_Zty(wheat_Zmat, y)
    x = gauss_seidel(A, b)
    assert isinstance(x, numpy.ndarray)
    assert x.shape == (wheat_nmkr,)

def test_rrBLUP_ML0(wheat_yvec, wheat_Zmat):
    out = rrBLUP_ML0(wheat_yvec, wheat_Zmat)
    assert isinstance(out, dict)

############################## Test class docstring ############################
def test_class_docstring():
    assert_docstring(rrBLUPModel0)

############################ Test Class Properties #############################

def test_u_misc_fget(algmod, algmod_u_misc):
    assert numpy.all(algmod.u_misc == algmod_u_misc)

def test_u_a_fget(algmod, algmod_u_a):
    assert numpy.all(algmod.u_a == algmod_u_a)

def test_beta_fget(algmod, algmod_beta):
    assert numpy.all(algmod.beta == algmod_beta)

def test_trait_fget(algmod, algmod_trait):
    assert numpy.all(algmod.trait == algmod_trait)

def test_algmod_model_name_fget(algmod, algmod_model_name):
    assert algmod.model_name == algmod_model_name

def test_algmod_params_fget(algmod, algmod_params):
    assert algmod.hyperparams == algmod_params

########################## Test Class Special Methods ##########################

### __init__
def test___init___is_concrete():
    assert_concrete_method(rrBLUPModel0, "__init__")

### __copy__
def test___copy___is_concrete():
    assert_concrete_method(rrBLUPModel0, "__copy__")

def test___copy__(algmod):
    tmp = algmod.__copy__()
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams
    tmp = copy.copy(algmod)
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams

### __deepcopy__
def test___deepcopy___is_concrete():
    assert_concrete_method(rrBLUPModel0, "__deepcopy__")

def test___deepcopy__(algmod):
    tmp = algmod.__deepcopy__()
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams
    tmp = copy.deepcopy(algmod)
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams

############################# Test concrete methods ############################

############# Copy methods #############

### copy
def test_copy_is_concrete():
    assert_concrete_method(rrBLUPModel0, "copy")

def test_copy(algmod):
    tmp = algmod.copy()
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams

### deepcopy
def test_deepcopy_is_concrete():
    assert_concrete_method(rrBLUPModel0, "deepcopy")

def test_deepcopy(algmod):
    tmp = algmod.deepcopy()
    assert numpy.all(tmp.beta == algmod.beta)
    assert numpy.all(tmp.u == algmod.u)
    assert numpy.all(tmp.u_misc == algmod.u_misc)
    assert numpy.all(tmp.u_a == algmod.u_a)
    assert numpy.all(tmp.trait == algmod.trait)
    assert tmp.model_name == algmod.model_name
    assert tmp.hyperparams == algmod.hyperparams

########## Prediction methods ##########

### fit_numpy
def test_fit_numpy_is_concrete():
    assert_concrete_method(rrBLUPModel0, "fit_numpy")

def test_fit_numpy(wheat_Ymat, wheat_Zmat):
    model = rrBLUPModel0.fit_numpy(wheat_Ymat, None, wheat_Zmat)
    assert isinstance(model, rrBLUPModel0)

### fit
def test_fit_is_concrete():
    assert_concrete_method(rrBLUPModel0, "fit")

def test_fit(bvmat, pgmat):
    model = rrBLUPModel0.fit(bvmat, None, pgmat)
    assert isinstance(model, rrBLUPModel0)

######################### Test class utility functions #########################

### check_is_rrBLUPModel0
def test_check_is_rrBLUPModel0_is_concrete():
    assert_concrete_function(check_is_rrBLUPModel0)

def test_check_is_rrBLUPModel0(algmod):
    with not_raises(TypeError):
        check_is_rrBLUPModel0(algmod, "algmod")
    with pytest.raises(TypeError):
        check_is_rrBLUPModel0(None, "algmod")
