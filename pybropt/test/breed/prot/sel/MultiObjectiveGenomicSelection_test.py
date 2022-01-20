import numpy
import pytest

from numpy.random import Generator
from numpy.random import PCG64

from pybropt.breed.prot.sel.MultiObjectiveGenomicSelection import MultiObjectiveGenomicSelection
from pybropt.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybropt.popgen.bvmat.DenseEstimatedBreedingValueMatrix import DenseEstimatedBreedingValueMatrix
from pybropt.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybropt.algo.opt.SteepestAscentSetHillClimber import SteepestAscentSetHillClimber
from pybropt.algo.opt.StochasticAscentSetHillClimber import StochasticAscentSetHillClimber
from pybropt.breed.prot.sel.transfn import trans_sum
from pybropt.breed.prot.sel.transfn import trans_ndpt_to_vec_dist

################################################################################
################################## Genotypes ###################################
################################################################################
@pytest.fixture
def mat_int8():
    yield numpy.int8([
       [[1, 0, 0, 0, 0, 0, 1, 0, 1, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 1]],
       [[0, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 1, 1, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0]]
    ])

@pytest.fixture
def mat_int8_big():
    yield numpy.int8([
       [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]],

       [[1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
        [1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]]
    ])

@pytest.fixture
def mat_chrgrp():
    yield numpy.int64([1, 1, 2, 2, 3, 3, 4, 4, 5, 5])

@pytest.fixture
def mat_chrgrp_big():
    yield numpy.int64([1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2])

@pytest.fixture
def mat_phypos():
    yield numpy.int64([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

@pytest.fixture
def mat_phypos_big():
    yield numpy.arange(16)

@pytest.fixture
def mat_taxa():
    yield numpy.object_(["Line1", "Line2", "Line3", "Line4", "Line5"])

@pytest.fixture
def mat_taxa_big():
    yield numpy.object_(["Line"+str(i).zfill(2) for i in range(30)])

@pytest.fixture
def mat_taxa_grp():
    yield numpy.int64([1, 1, 2, 2, 2])

@pytest.fixture
def dpgvmat(mat_int8, mat_chrgrp, mat_phypos, mat_taxa, mat_taxa_grp):
    yield DensePhasedGenotypeMatrix(
        mat = mat_int8,
        vrnt_chrgrp = mat_chrgrp,
        vrnt_phypos = mat_phypos,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp
    )

@pytest.fixture
def dpgvmat_big(mat_int8_big, mat_chrgrp_big, mat_phypos_big, mat_taxa_big):
    yield DensePhasedGenotypeMatrix(
        mat = mat_int8_big,
        vrnt_chrgrp = mat_chrgrp_big,
        vrnt_phypos = mat_phypos_big,
        taxa = mat_taxa_big
    )

################################################################################
################################ Genomic model #################################
################################################################################
@pytest.fixture
def beta():
    yield numpy.float64([
        [1.4],
        [2.5],
        [7.2]
    ])

@pytest.fixture
def u():
    yield numpy.float64([
        [-0.33,  2.08, -2.42],
        [-0.69, -1.87, -1.38],
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
def u_big():
    yield numpy.float64([
       [-0.87, -0.16, -0.04],
       [-0.03,  0.05, -0.15],
       [ 0.36, -0.15,  0.54],
       [ 2.35,  1.12,  0.33],
       [-0.93, -0.59, -1.1 ],
       [-0.63,  0.61,  0.15],
       [-0.8 ,  0.95, -0.56],
       [-1.03,  1.55,  0.12],
       [-1.21,  0.79,  1.42],
       [-0.09, -1.88, -1.83],
       [-0.03,  1.97, -1.98],
       [-0.04,  0.2 ,  1.43],
       [ 0.45,  1.26, -2.21],
       [-0.31, -0.62,  1.09],
       [ 0.9 , -1.37,  0.91],
       [-0.87,  1.2 , -1.68]
   ])

@pytest.fixture
def trait():
    yield numpy.object_(["protein", "yield", "quality"])

@pytest.fixture
def model_name():
    yield "test_glgmod"

@pytest.fixture
def params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def glgmod(beta, u, trait, model_name, params):
    yield AdditiveLinearGenomicModel(
        beta = beta,
        u = u,
        trait = trait,
        model_name = model_name,
        params = params
    )

@pytest.fixture
def glgmod_big(beta, u_big, trait, model_name, params):
    yield AdditiveLinearGenomicModel(
        beta = beta,
        u = u_big,
        trait = trait,
        model_name = model_name,
        params = params
    )

################################################################################
############################ Breeding values model #############################
################################################################################
@pytest.fixture
def bvmat(glgmod, dpgvmat):
    yield glgmod.predict(dpgvmat)

@pytest.fixture
def bvmat_big(glgmod_big, dpgvmat_big):
    yield glgmod_big.predict(dpgvmat_big)

################################################################################
##################### MultiObjectiveGenomicSelection #####################
################################################################################
@pytest.fixture
def nparent():
    yield 2

@pytest.fixture
def ncross():
    yield 1

@pytest.fixture
def nprogeny():
    yield 10

@pytest.fixture
def method():
    yield "pareto"

@pytest.fixture
def objfn_trans():
    yield trans_sum

@pytest.fixture
def objfn_wt():
    yield numpy.array([-1.0,-1.0])  # all objectives are minimizing

@pytest.fixture
def ndset_trans():
    yield trans_ndpt_to_vec_dist

@pytest.fixture
def ndset_trans_kwargs(objfn_wt):
    yield {
        "objfn_wt" : numpy.array([0.3,0.7]),
        "wt" : objfn_wt
    }

@pytest.fixture
def ndset_wt():
    yield -1.0    # function is minimizing (minimize distance)

@pytest.fixture
def rng():
    yield Generator(PCG64(192837465))

@pytest.fixture
def algorithm(nparent, dpgvmat, rng):
    yield StochasticAscentSetHillClimber(
        k = nparent,
        setspace = numpy.arange(dpgvmat.ntaxa),
        rng = rng,
        objwt = 1.0,
    )

@pytest.fixture
def mogps(nparent, ncross, nprogeny, algorithm, method, objfn_trans, objfn_wt, ndset_trans, ndset_trans_kwargs, ndset_wt, rng):
    yield MultiObjectiveGenomicSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        algorithm = algorithm,
        method = method,
        objfn_trans = objfn_trans,
        objfn_wt = objfn_wt,
        ndset_trans = ndset_trans,
        ndset_trans_kwargs = ndset_trans_kwargs,
        ndset_wt = ndset_wt,
        rng = rng
    )

################################################################################
#################################### Tests #####################################
################################################################################

########################################
########### test calc_mkrwt ############
########################################
def test_calc_mkrwt_magnitude(u_big):
    a = MultiObjectiveGenomicSelection._calc_mkrwt("magnitude", u_big)
    b = numpy.absolute(u_big)
    assert numpy.all(a == b)

def test_calc_mkrwt_equal(u_big):
    a = MultiObjectiveGenomicSelection._calc_mkrwt("equal", u_big)
    assert numpy.all(a == 1.0)

def test_calc_mkrwt_str_case(u_big):
    a = MultiObjectiveGenomicSelection._calc_mkrwt("mAgNiTuDe", u_big)
    b = numpy.absolute(u_big)
    assert numpy.all(a == b)
    a = MultiObjectiveGenomicSelection._calc_mkrwt("Equal", u_big)
    assert numpy.all(a == 1.0)

def test_calc_mkrwt_str_ValueError(u_big):
    with pytest.raises(ValueError):
        a = MultiObjectiveGenomicSelection._calc_mkrwt("unknown", u_big)

def test_calc_mkrwt_ndarray(u_big):
    wt = numpy.random.normal(size = u_big.shape)
    a = MultiObjectiveGenomicSelection._calc_mkrwt(wt, u_big)
    assert numpy.all(a == wt)

def test_calc_mkrwt_type_TypeError(u_big):
    with pytest.raises(TypeError):
        a = MultiObjectiveGenomicSelection._calc_mkrwt(None, u_big)

########################################
########### test calc_tfreq ############
########################################
def test_calc_tfreq_positive(u_big):
    a = MultiObjectiveGenomicSelection._calc_tfreq("positive", u_big)
    b = numpy.float64(u_big >= 0.0)
    assert numpy.all(a == b)

def test_calc_tfreq_negative(u_big):
    a = MultiObjectiveGenomicSelection._calc_tfreq("negative", u_big)
    b = numpy.float64(u_big <= 0.0)
    assert numpy.all(a == b)

def test_calc_tfreq_stabilizing(u_big):
    a = MultiObjectiveGenomicSelection._calc_tfreq("stabilizing", u_big)
    assert numpy.all(a == 0.5)

def test_calc_tfreq_str_case(u_big):
    a = MultiObjectiveGenomicSelection._calc_tfreq("PoSiTiVe", u_big)
    b = numpy.float64(u_big >= 0.0)
    assert numpy.all(a == b)
    a = MultiObjectiveGenomicSelection._calc_tfreq("NEGATIVE", u_big)
    b = numpy.float64(u_big <= 0.0)
    assert numpy.all(a == b)
    a = MultiObjectiveGenomicSelection._calc_tfreq("Stabilizing", u_big)
    assert numpy.all(a == 0.5)

def test_calc_tfreq_str_ValueError(u_big):
    with pytest.raises(ValueError):
        a = MultiObjectiveGenomicSelection._calc_tfreq("unknown", u_big)

def test_calc_tfreq_ndarray(u_big):
    wt = numpy.random.uniform(0, 1, size = u_big.shape)
    a = MultiObjectiveGenomicSelection._calc_tfreq(wt, u_big)
    assert numpy.all(a == wt)

def test_calc_tfreq_type_TypeError(u_big):
    with pytest.raises(TypeError):
        a = MultiObjectiveGenomicSelection._calc_tfreq(None, u_big)

# test constructor
def test_init(mogps):
    assert True

# test single objective optimization
def test_pselect_single(mogps, dpgvmat, bvmat, glgmod, ncross, nprogeny):
    geno = {
        "cand" : dpgvmat,
        "main" : dpgvmat,
        "queue" : [dpgvmat]
    }
    bval = {
        "cand" : bvmat,
        "cand_true" : bvmat,
        "main" : bvmat,
        "main_true" : bvmat
    }
    gmod = {
        "cand" : glgmod,
        "main" : glgmod,
        "true" : glgmod
    }

    out_gmat, out_sel, out_ncross, out_nprogeny, out_misc = mogps.pselect(
        t_cur = 0,
        t_max = 20,
        geno = geno,
        bval = bval,
        gmod = gmod,
        method = "single"
    )

    # should have objfn_eval = [0., 24.085] = [perfect score, min dist]
    assert numpy.all(out_sel == [0,3]) or numpy.all(out_sel == [3,0])
    assert out_ncross == ncross
    assert out_nprogeny == nprogeny

def test_ppareto(mogps, dpgvmat_big, bvmat_big, glgmod_big, ncross, nprogeny):
    geno = {
        "cand" : dpgvmat_big,
        "main" : dpgvmat_big,
        "queue" : [dpgvmat_big]
    }
    bval = {
        "cand" : bvmat_big,
        "cand_true" : bvmat_big,
        "main" : bvmat_big,
        "main_true" : bvmat_big
    }
    gmod = {
        "cand" : glgmod_big,
        "main" : glgmod_big,
        "true" : glgmod_big
    }

    frontier, pop, logbook = mogps.ppareto(
        t_cur = 0,
        t_max = 20,
        geno = geno,
        bval = bval,
        gmod = gmod,
        nparent = 6
    )

    assert isinstance(frontier, numpy.ndarray)

    from matplotlib import pyplot

    pyplot.scatter(frontier[:,0], frontier[:,1], c="b")
    pyplot.axis("tight")
    pyplot.savefig("frontier_new.png")

def test_pselect_pareto(mogps, dpgvmat_big, bvmat_big, glgmod_big, ncross, nprogeny):
    geno = {
        "cand" : dpgvmat_big,
        "main" : dpgvmat_big,
        "queue" : [dpgvmat_big]
    }
    bval = {
        "cand" : bvmat_big,
        "cand_true" : bvmat_big,
        "main" : bvmat_big,
        "main_true" : bvmat_big
    }
    gmod = {
        "cand" : glgmod_big,
        "main" : glgmod_big,
        "true" : glgmod_big
    }

    out_gmat, out_sel, out_ncross, out_nprogeny, out_misc = mogps.pselect(
        t_cur = 0,
        t_max = 20,
        geno = geno,
        bval = bval,
        gmod = gmod,
        nparent = 6,
        method = "pareto"
    )

    # true solution = [28 17  0  8 23 29]
    soln = [28, 17, 0, 8, 23, 29]
    assert numpy.all(numpy.in1d(soln, out_sel))
    assert out_ncross == ncross
    assert out_nprogeny == nprogeny
