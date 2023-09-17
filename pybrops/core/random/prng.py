"""
Module containing wrappers for random number generation.
"""

# list of all names to import when using wildcard ('*') import
__all__ = [
    "global_prng",
    "beta",
    "binomial",
    "bytes",
    "chisquare",
    "choice",
    "dirichlet",
    "exponential",
    "f",
    "gamma",
    "geometric",
    "gumbel",
    "hypergeometric",
    "laplace",
    "logistic",
    "lognormal",
    "logseries",
    "multinomial",
    "multivariate_normal",
    "negative_binomial",
    "noncentral_chisquare",
    "noncentral_f",
    "normal",
    "pareto",
    "permutation",
    "poisson",
    "power",
    # "rand",
    # "randint",
    # "randn",
    "random",
    # "random_integers",
    # "random_sample",
    "rayleigh",
    "shuffle",
    "spawn",
    "standard_cauchy",
    "standard_exponential",
    "standard_gamma",
    "standard_normal",
    "standard_t",
    "triangular",
    "uniform",
    "vonmises",
    "wald",
    "weibull",
    "zipf",
    # "sample",
    # "ranf",
    "seed",
    "spawn",
]

from typing import Optional, Union
import numpy
import random as py_random
from numpy.random import PCG64
from numpy.random import Generator

# get the numpy random number generator state
global_prng = numpy.random.random.__self__

# expose methods of global_prng to the public
beta = global_prng.beta
binomial = global_prng.binomial
bytes = global_prng.bytes
chisquare = global_prng.chisquare
choice = global_prng.choice
dirichlet = global_prng.dirichlet
exponential = global_prng.exponential
f = global_prng.f
gamma = global_prng.gamma
geometric = global_prng.geometric
gumbel = global_prng.gumbel
hypergeometric = global_prng.hypergeometric
laplace = global_prng.laplace
logistic = global_prng.logistic
lognormal = global_prng.lognormal
logseries = global_prng.logseries
multinomial = global_prng.multinomial
multivariate_normal = global_prng.multivariate_normal
negative_binomial = global_prng.negative_binomial
noncentral_chisquare = global_prng.noncentral_chisquare
noncentral_f = global_prng.noncentral_f
normal = global_prng.normal
pareto = global_prng.pareto
permutation = global_prng.permutation
poisson = global_prng.poisson
power = global_prng.power
# rand = global_prng.rand
# randint = global_prng.randint
# randn = global_prng.randn
random = global_prng.random
# random_integers = global_prng.random_integers
# random_sample = global_prng.random_sample
rayleigh = global_prng.rayleigh
shuffle = global_prng.shuffle
standard_cauchy = global_prng.standard_cauchy
standard_exponential = global_prng.standard_exponential
standard_gamma = global_prng.standard_gamma
standard_normal = global_prng.standard_normal
standard_t = global_prng.standard_t
triangular = global_prng.triangular
uniform = global_prng.uniform
vonmises = global_prng.vonmises
wald = global_prng.wald
weibull = global_prng.weibull
zipf = global_prng.zipf
# Two legacy that are trivial wrappers around random_sample
# sample = global_prng.random_sample
# ranf = global_prng.random_sample

# random number seeding
def seed(s: Optional[int] = None) -> None:
    """
    Seed the global pseudo-random number generator (PRNG). Begins by seeding the
    ``random`` module, then seeds ``numpy.random``, which combined finishes the
    seeding process for ``pybrops.core.random``. This ensures that both the
    ``random`` and ``numpy.random`` generators are tied to the same entropy
    source.

    Parameters
    ----------
    s : None, int
        Seed to use for the PRNG. If None, use either system time, or random
        bits from ``os.urandom``. This behavior depends entirely on the
        system implementation of ``random.seed()``.
    """
    py_random.seed(s)                                   # seed random module
    numpy.random.seed(py_random.randint(0, 2**32-1))    # seed numpy.random with 4 bytes of entropy

# random number generator spawner
def spawn(
        n: Optional[int] = None, 
        BitGenerator: numpy.random.BitGenerator = PCG64, 
        sbits: int = 64
    ) -> Union[numpy.random.Generator,list]:
    """
    Spawn new numpy PRNG streams.

    Parameters
    ----------
    n : None, int
        Number of streams to spawn using entropy from the ``random`` module.
        If ``None``, generate one stream. Integer must be positive or zero.
    BitGenerator : numpy.random.BitGenerator
        A constructor for a ``numpy.random.BitGenerator`` object.
    sbits : int
        Number of bits used to seed the new PRNG streams.

    Returns
    -------
    streams : numpy.random.Generator, list of numpy.random.Generator
        New random number generator streams.
    """
    if n is None:
        # seed a single generator using 'sbits' of entropy
        out = Generator(BitGenerator(py_random.randint(0, 2**sbits-1)))
    elif isinstance(n, int):
        if n < 0:
            raise ValueError("'{0}' must be greater than or equal to 0: {0} = {1}".format("n", n))
        # seed 'n' generators using 'sbits' of entropy
        out = [Generator(BitGenerator(py_random.randint(0, 2**sbits-1))) for _ in range(n)]
    else:
        raise TypeError("'{0}' must be of type int".format("n"))
    return out
