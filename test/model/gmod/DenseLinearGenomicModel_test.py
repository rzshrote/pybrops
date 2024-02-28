from typing import Dict, Optional
import numpy
from pandas.core.api import DataFrame as DataFrame
import pytest
import os.path
from pybrops.core.io.CSVDictInputOutput import CSVDictInputOutput
from pybrops.core.io.PandasDictInputOutput import PandasDictInputOutput
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

from pybrops.test.assert_python import assert_classmethod_isconcrete, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.gmod.DenseLinearGenomicModel import DenseLinearGenomicModel
from pybrops.model.gmod.DenseLinearGenomicModel import check_is_DenseLinearGenomicModel

from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix

################################################################################
################################ Test fixtures #################################
################################################################################

class DummyDenseLinearGenomicModel(DenseLinearGenomicModel):
    __doc__ = DenseLinearGenomicModel.__doc__
    def copy(self) -> GenomicModel:
        return super().copy()
    def deepcopy(self, memo: dict | None) -> GenomicModel:
        return super().deepcopy(memo)
    @classmethod
    def from_csv_dict(cls, filenames: Dict[str, str], **kwargs: dict) -> CSVDictInputOutput:
        return super().from_csv_dict(filenames, **kwargs)
    @classmethod
    def from_pandas_dict(cls, dic: Dict[str, DataFrame], **kwargs: dict) -> PandasDictInputOutput:
        return super().from_pandas_dict(dic, **kwargs)
    def gegv(self, gtobj: GenotypeMatrix, **kwargs: dict) -> BreedingValueMatrix:
        return super().gegv(gtobj, **kwargs)
    def gegv_numpy(self, Z: numpy.ndarray, **kwargs: dict) -> numpy.ndarray:
        return super().gegv_numpy(Z, **kwargs)
    @property
    def nexplan(self) -> object:
        return super().nexplan
    @property
    def nexplan_beta(self) -> object:
        return super().nexplan_beta
    @property
    def nexplan_u(self) -> object:
        return super().nexplan_u
    @property
    def nparam(self) -> object:
        return super().nparam
    def to_csv_dict(self, filenames: Dict[str, str], **kwargs: dict) -> None:
        return super().to_csv_dict(filenames, **kwargs)
    def to_pandas_dict(self, **kwargs: dict) -> Dict[str, DataFrame]:
        return super().to_pandas_dict(**kwargs)

############################################################
###################### Genomic model #######################
############################################################
@pytest.fixture
def mat_beta():
    yield numpy.float64([
        [1.4, 2.5, 7.2]
    ])

@pytest.fixture
def mat_u():
    yield numpy.float64([
        [-0.33,  2.08, -2.42],
        [-0.69, -1.87,  1.38],
        [ 1.12,  1.38, -5.65],
        [-1.44,  0.20,  4.22],
        [ 0.88, -0.81,  1.55],
        [ 1.23,  0.25,  5.13],
        [ 0.19,  4.35,  0.15],
        [-2.12,  0.73, -0.38],
        [-0.87,  1.25,  2.38],
        [ 0.06, -2.52,  2.48]
    ])

@pytest.fixture
def mat_trait():
    yield numpy.object_(["protein", "yield", "quality"])

@pytest.fixture
def model_name():
    yield "test_glgmod"

@pytest.fixture
def hyperparams():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def glgmod(mat_beta, mat_u, mat_trait, model_name, hyperparams):
    yield DummyDenseLinearGenomicModel(
        beta = mat_beta,
        u = mat_u,
        trait = mat_trait,
        model_name = model_name,
        hyperparams = hyperparams
    )

############################################################
######################## Genotypes #########################
############################################################
@pytest.fixture
def mat_int8():
    yield numpy.int8([
       [[1, 0, 1, 0, 0, 0, 1, 0, 1, 1],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 1, 1, 0, 0, 0, 1, 1],
        [1, 0, 1, 0, 0, 1, 0, 0, 0, 1]],
       [[0, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        [0, 0, 1, 0, 1, 0, 0, 1, 1, 0],
        [0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 1, 1, 0, 0, 0, 0, 0, 0]]
    ])

@pytest.fixture
def mat_chrgrp():
    yield numpy.int64([1, 1, 2, 2, 3, 3, 4, 4, 5, 5])

@pytest.fixture
def mat_phypos():
    yield numpy.int64([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

@pytest.fixture
def mat_taxa():
    yield numpy.object_(["Line1", "Line2", "Line3", "Line4", "Line5"])

@pytest.fixture
def mat_taxa_grp():
    yield numpy.int64([1, 1, 2, 2, 2])

@pytest.fixture
def dpgmat(mat_int8, mat_chrgrp, mat_phypos, mat_taxa, mat_taxa_grp):
    yield DensePhasedGenotypeMatrix(
        mat = mat_int8,
        vrnt_chrgrp = mat_chrgrp,
        vrnt_phypos = mat_phypos,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp
    )

############################################################
##################### Breeding values ######################
############################################################
@pytest.fixture
def mat_intercept(dpgmat, mat_beta):
    n = dpgmat.ntaxa
    q = mat_beta.shape[0]
    a = numpy.ones((n,1), dtype = "float64")
    yield a

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(DummyDenseLinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "__init__")

def test_copy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "__copy__")

def test_deepcopy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "__deepcopy__")

def test_fit_numpy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "fit_numpy")

def test_fit_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "fit")

def test_predict_numpy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "predict_numpy")

def test_predict_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "predict")

def test_score_numpy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "score_numpy")

def test_score_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "score")

def test_gebv_numpy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "gebv_numpy")

def test_gebv_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "gebv")

def test_var_G_numpy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "var_G_numpy")

def test_var_G_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "var_G")

def test_var_A_numpy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "var_A_numpy")

def test_var_A_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "var_A")

def test_var_a_numpy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "var_a_numpy")

def test_var_a_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "var_a")

def test_bulmer_numpy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "bulmer_numpy")

def test_bulmer_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "bulmer")

def test_usl_numpy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "usl_numpy")

def test_usl_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "usl")

def test_lsl_numpy_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "lsl_numpy")

def test_lsl_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "lsl")

def test_from_hdf5_is_concrete():
    assert_classmethod_isconcrete(DummyDenseLinearGenomicModel, "from_hdf5")

def test_to_hdf5_is_concrete():
    assert_method_isconcrete(DummyDenseLinearGenomicModel, "to_hdf5")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_u_fget(glgmod, mat_u):
    assert numpy.all(glgmod.u == mat_u)

def test_beta_fget(glgmod, mat_beta):
    assert numpy.all(glgmod.beta == mat_beta)

def test_trait_fget(glgmod, mat_trait):
    assert numpy.all(glgmod.trait == mat_trait)

def test_model_name_fget(glgmod, model_name):
    assert glgmod.model_name == model_name

def test_params_fget(glgmod, hyperparams):
    assert glgmod.hyperparams == hyperparams

################################################################################
###################### Test concrete method functionality ######################
################################################################################

########################################
########### Prediction tests ###########
########################################
def test_fit_numpy(glgmod):
    with pytest.raises(AttributeError):
        glgmod.fit_numpy(None, None, None)

def test_fit(glgmod):
    with pytest.raises(AttributeError):
        glgmod.fit(None, None, None)

def test_predict_numpy(glgmod, mat_intercept, mat_beta, mat_int8, mat_u):
    geno = mat_int8.sum(0)
    a = glgmod.predict_numpy(mat_intercept, geno)
    b = (mat_intercept @ mat_beta) + (geno @ mat_u)
    assert numpy.all(a == b)
    assert isinstance(a, numpy.ndarray)

def test_predict(glgmod, mat_intercept, mat_beta, mat_int8, mat_u, dpgmat):
    geno = mat_int8.sum(0)
    a = glgmod.predict(mat_intercept, dpgmat)
    b = (mat_intercept @ mat_beta) + (geno @ mat_u)
    b = (b - b.mean(0)) / b.std(0)
    assert numpy.all(a == b)
    assert isinstance(a, BreedingValueMatrix)

def test_score_numpy(glgmod, mat_intercept, mat_int8):
    geno = mat_int8.sum(0)
    y_true = glgmod.predict_numpy(mat_intercept, geno)
    out = glgmod.score_numpy(y_true, mat_intercept, geno)
    assert isinstance(out, numpy.ndarray)
    assert numpy.all(out == 1.0)
    assert len(out) == glgmod.ntrait

def test_score(glgmod, mat_intercept, dpgmat):
    y_true = glgmod.predict(mat_intercept, dpgmat)
    out = glgmod.score(y_true, mat_intercept, dpgmat)
    assert isinstance(out, numpy.ndarray)
    assert numpy.all(out == 1.0)

########################################
###### Variance calculation tests ######
########################################


########################################
######## Selection limit tests #########
########################################
def test_usl(glgmod, dpgmat, mat_u, mat_int8):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = mat_u.T[:,None,None,:] * mat_int8
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.max((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)
    a_usl = uhaplo.sum(1)
    b_usl = glgmod.usl(dpgmat)

    assert numpy.all(a_usl == b_usl)

def test_usl(glgmod, dpgmat, mat_u, mat_int8):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = mat_u.T[:,None,None,:] * mat_int8
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.min((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)    pass

    a_lsl = uhaplo.sum(1)
    b_lsl = glgmod.lsl(dpgmat)

    assert numpy.all(a_lsl == b_lsl)

### File I/O tests ###
def test_to_from_hdf5(glgmod, shared_datadir):
    glgmod.to_hdf5(shared_datadir / "glgmod.hdf5")
    glgmod.to_hdf5(shared_datadir / "glgmod.hdf5", "prefix")

    # test whether file was created
    assert os.path.isfile(shared_datadir / "glgmod.hdf5")

    glgmod1 = DummyDenseLinearGenomicModel.from_hdf5(shared_datadir / "glgmod.hdf5")
    glgmod2 = DummyDenseLinearGenomicModel.from_hdf5(
        shared_datadir / "glgmod.hdf5",
        "prefix"
    )

    # test whether data was loaded properly
    assert numpy.all(glgmod.beta == glgmod1.beta)
    assert numpy.all(glgmod.u == glgmod1.u)
    assert numpy.all(glgmod.trait == glgmod1.trait)
    assert glgmod.model_name == glgmod1.model_name
    assert glgmod.hyperparams == glgmod1.hyperparams

    assert numpy.all(glgmod.beta == glgmod2.beta)
    assert numpy.all(glgmod.u == glgmod2.u)
    assert numpy.all(glgmod.trait == glgmod2.trait)
    assert glgmod.model_name == glgmod2.model_name
    assert glgmod.hyperparams == glgmod2.hyperparams

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_DenseLinearGenomicModel_is_concrete():
    assert_function_isconcrete(check_is_DenseLinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseLinearGenomicModel(glgmod):
    with not_raises(TypeError):
        check_is_DenseLinearGenomicModel(glgmod, "glgmod")
    with pytest.raises(TypeError):
        check_is_DenseLinearGenomicModel(None, "glgmod")
