import inspect
import pytest
import numpy
from numpy.random import Generator

from pybrops.core.random import (
    default_BitGenerator,
    seed_nbytes,
    beta,
    binomial,
    bytes,
    chisquare,
    choice,
    dirichlet,
    exponential,
    f,
    gamma,
    geometric,
    gumbel,
    hypergeometric,
    integers,
    laplace,
    logistic,
    lognormal,
    logseries,
    multinomial,
    multivariate_hypergeometric,
    multivariate_normal,
    negative_binomial,
    noncentral_chisquare,
    noncentral_f,
    normal,
    pareto,
    permutation,
    poisson,
    power,
    random,
    rayleigh,
    seed,
    shuffle,
    spawn,
    standard_cauchy,
    standard_exponential,
    standard_gamma,
    standard_normal,
    standard_t,
    triangular,
    uniform,
    vonmises,
    wald,
    weibull,
    zipf,
)

################################################################################
############################# Assert import usage ##############################
################################################################################

def test_pybrops_core_random_import():
    import pybrops.core.random
    rng = None

    if rng is None:
        rng = pybrops.core.random

    s = (10,10)
    m = rng.beta(3.1,5.2, s)
    assert m.shape == s

################################################################################
##################### Test that all functions are callable #####################
################################################################################

def test_default_BitGenerator_callable():
    assert callable(default_BitGenerator)

def test_beta_callable():
    assert callable(beta)

def test_binomial_callable():
    assert callable(binomial)

def test_bytes_callable():
    assert callable(bytes)

def test_chisquare_callable():
    assert callable(chisquare)

def test_choice_callable():
    assert callable(choice)

def test_dirichlet_callable():
    assert callable(dirichlet)

def test_exponential_callable():
    assert callable(exponential)

def test_f_callable():
    assert callable(f)

def test_gamma_callable():
    assert callable(gamma)

def test_geometric_callable():
    assert callable(geometric)

def test_gumbel_callable():
    assert callable(gumbel)

def test_hypergeometric_callable():
    assert callable(hypergeometric)

def test_integers_callable():
    assert callable(integers)

def test_laplace_callable():
    assert callable(laplace)

def test_logistic_callable():
    assert callable(logistic)

def test_lognormal_callable():
    assert callable(lognormal)

def test_logseries_callable():
    assert callable(logseries)

def test_multinomial_callable():
    assert callable(multinomial)

def test_multivariate_hypergeometric_callable():
    assert callable(multivariate_hypergeometric)

def test_multivariate_normal_callable():
    assert callable(multivariate_normal)

def test_negative_binomial_callable():
    assert callable(negative_binomial)

def test_noncentral_chisquare_callable():
    assert callable(noncentral_chisquare)

def test_noncentral_f_callable():
    assert callable(noncentral_f)

def test_normal_callable():
    assert callable(normal)

def test_pareto_callable():
    assert callable(pareto)

def test_permutation_callable():
    assert callable(permutation)

def test_poisson_callable():
    assert callable(poisson)

def test_power_callable():
    assert callable(power)

def test_random_callable():
    assert callable(random)

def test_rayleigh_callable():
    assert callable(rayleigh)

def test_seed_callable():
    assert callable(seed)

def test_shuffle_callable():
    assert callable(shuffle)

def test_spawn_callable():
    assert callable(spawn)

def test_standard_cauchy_callable():
    assert callable(standard_cauchy)

def test_standard_exponential_callable():
    assert callable(standard_exponential)

def test_standard_gamma_callable():
    assert callable(standard_gamma)

def test_standard_normal_callable():
    assert callable(standard_normal)

def test_standard_t_callable():
    assert callable(standard_t)

def test_triangular_callable():
    assert callable(triangular)

def test_uniform_callable():
    assert callable(uniform)

def test_vonmises_callable():
    assert callable(vonmises)

def test_wald_callable():
    assert callable(wald)

def test_weibull_callable():
    assert callable(weibull)

def test_zipf_callable():
    assert callable(zipf)

################################################################################
############################ Test seeding/spawning #############################
################################################################################

def test_spawn_None():
    s = spawn()
    assert isinstance(s, Generator)

def test_spawn_int():
    l = spawn(5)
    assert isinstance(l, list)
    for e in l:
        assert isinstance(e, Generator)
    s = (10,10)
    m = beta(3.1, 5.2, s)
    assert m.shape == s

def test_seed_None():
    seed(None)

def test_seed_int():
    seed(int(123456789))

################################################################################
######################## Test random number generation #########################
################################################################################

def test_beta_random():
    s = (10,10)
    m = beta(1, 1, size = s)
    assert m.shape == s

def test_binomial_random():
    s = (10,10)
    m = binomial(10, 0.5, size = s)
    assert m.shape == s

def test_bytes_random():
    s = 10
    m = bytes(s)
    assert len(m) == s

def test_chisquare_random():
    s = (10,10)
    m = chisquare(1, size = s)
    assert m.shape == s

def test_choice_random():
    s = (10,10)
    m = choice([0,1,2,3,4], size = s, replace = True)
    assert m.shape == s

def test_dirichlet_random():
    s = (10,10)
    m = dirichlet([0.5,0.5], size = s)
    assert m.shape == (s + (2,))

def test_exponential_random():
    s = (10,10)
    m = exponential(1, size = s)
    assert m.shape == s

def test_f_random():
    s = (10,10)
    m = f(1, 1, size = s)
    assert m.shape == s

def test_gamma_random():
    s = (10,10)
    m = gamma(1, 1, size = s)
    assert m.shape == s

def test_geometric_random():
    s = (10,10)
    m = geometric(0.5, size = s)
    assert m.shape == s

def test_gumbel_random():
    s = (10,10)
    m = gumbel(1, 1, size = s)
    assert m.shape == s

def test_hypergeometric_random():
    s = (10,10)
    m = hypergeometric(10, 10, 10, size = s)
    assert m.shape == s

# TODO: implement test
# def test_integers_random():
#     s = (10,10)
#     m = integers(1, 1, size = s)
#     assert m.shape == s

def test_laplace_random():
    s = (10,10)
    m = laplace(1, 1, size = s)
    assert m.shape == s

def test_logistic_random():
    s = (10,10)
    m = logistic(1, 1, size = s)
    assert m.shape == s

def test_lognormal_random():
    s = (10,10)
    m = lognormal(1, 1, size = s)
    assert m.shape == s

def test_logseries_random():
    s = (10,10)
    m = logseries(0.5, size = s)
    assert m.shape == s

def test_multinomial_random():
    s = (10,10)
    m = multinomial(10, [0.25,0.25,0.25,0.25], size = s)
    assert m.shape == (s + (4,))

# TODO: implement test
# def test_multivariate_hypergeometric_random():
#     s = (10,10)
#     m = multivariate_hypergeometric(1, 1, size = s)
#     assert m.shape == s

# TODO: implement test
# def test_multivariate_normal_random():
#     s = (10,10)
#     m = multivariate_normal(1, 1, size = s)
#     assert m.shape == s

def test_negative_binomial_random():
    s = (10,10)
    m = negative_binomial(10, 0.5, size = s)
    assert m.shape == s

def test_noncentral_chisquare_random():
    s = (10,10)
    m = noncentral_chisquare(1, 1, size = s)
    assert m.shape == s

def test_noncentral_f_random():
    s = (10,10)
    m = noncentral_f(1, 1, 1, size = s)
    assert m.shape == s

def test_normal_random():
    s = (10,10)
    m = normal(1, 1, size = s)
    assert m.shape == s

def test_pareto_random():
    s = (10,10)
    m = pareto(1, size = s)
    assert m.shape == s

def test_permutation_random():
    s = 10
    m = permutation(s)
    assert len(m) == s

def test_poisson_random():
    s = (10,10)
    m = poisson(1, size = s)
    assert m.shape == s

def test_power_random():
    s = (10,10)
    m = power(1, size = s)
    assert m.shape == s

def test_random_random():
    s = (10,10)
    m = random(size = s)
    assert m.shape == s

def test_rayleigh_random():
    s = (10,10)
    m = rayleigh(1, size = s)
    assert m.shape == s

def test_shuffle_random():
    s = 10
    m = numpy.arange(s)
    shuffle(m)
    assert len(m) == s

def test_standard_cauchy_random():
    s = (10,10)
    m = standard_cauchy(size = s)
    assert m.shape == s

def test_standard_exponential_random():
    s = (10,10)
    m = standard_exponential(size = s)
    assert m.shape == s

def test_standard_gamma_random():
    s = (10,10)
    m = standard_gamma(1, size = s)
    assert m.shape == s

def test_standard_normal_random():
    s = (10,10)
    m = standard_normal(size = s)
    assert m.shape == s

def test_standard_t_random():
    s = (10,10)
    m = standard_t(1, size = s)
    assert m.shape == s

def test_triangular_random():
    s = (10,10)
    m = triangular(1, 1.5, 2, size = s)
    assert m.shape == s

def test_uniform_random():
    s = (10,10)
    m = uniform(1, 10, size = s)
    assert m.shape == s

def test_vonmises_random():
    s = (10,10)
    m = vonmises(1, 1, size = s)
    assert m.shape == s

def test_wald_random():
    s = (10,10)
    m = wald(1, 1, size = s)
    assert m.shape == s

def test_weibull_random():
    s = (10,10)
    m = weibull(1, size = s)
    assert m.shape == s

def test_zipf_random():
    s = (10,10)
    m = zipf(10, size = s)
    assert m.shape == s
