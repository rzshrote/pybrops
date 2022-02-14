import random as py_random
import numpy
from numpy.random import PCG64
from numpy.random import Generator

class RandomStateWrapper:
    """
    Class to wrap seeding and random number generation functions.
    """
    def __init__(self, s = None):
        self.default_BitGenerator = PCG64   # default BitGenerator to use for random number generation
        self.seed_nbytes = 8                # number of bytes to use when seeding prng
        self.prng = None                    # actual generator object to use for creating random numbers
        self.seed(s)                        # seed (makes a Generator for prng)

    def spawn(self, n = None, BitGenerator = None):
        """
        Spawn new numpy PRNG streams.

        Parameters
        ----------
        n : None, int
            Number of streams to spawn using entropy from the random module.
            If None, generate one stream. Integer must be positive or zero.
        BitGenerator : function
            A function generating a BitGenerator object.

        Returns
        -------
        streams : numpy.random.Generator, list of numpy.random.Generator
            New random number generator streams.
        """
        # if no default bit generator constructor was provided, use default
        if BitGenerator is None:
            BitGenerator = self.default_BitGenerator

        if n is None:
            b = py_random.randbytes(self.seed_nbytes)   # generate random bytes
            s = int.from_bytes(b, byteorder = "big")    # convert bytes to int
            out = Generator(BitGenerator(s))            # create PRNG
        elif isinstance(n, int):
            if n < 0:
                raise ValueError("'{0}' must be greater than or equal to 0: {0} = {1}".format("n", n))
            bb = [py_random.randbytes(self.seed_nbytes) for _ in range(n)]
            ss = [int.from_bytes(b, byteorder = "big") for b in bb]
            out = [Generator(BitGenerator(s)) for s in ss]
        else:
            raise TypeError("'{0}' must be of type int".format("n"))
        return out

    def seed(self, s = None):
        """
        Seed a new pseudo-random number generator using the default BitGenerator.
        Begins by seeding the random module, then seeds numpy.random, then seeds
        pybrops.core.random itself.

        Parameters
        ----------
        s : None, int
            Seed to use for the PRNG. If None, use either system time, or random
            bits from os.urandom. This behavior depends entirely on the
            functionality of random.seed().
        """
        py_random.seed(s)                                   # seed random module
        s_byte = py_random.randbytes(4)                     # get 4 bytes of entropy
        s_int = int.from_bytes(s_byte, byteorder = 'big')   # convert entropy to integer
        numpy.random.seed(s_int)                            # seed numpy.random with 4 bytes of entropy
        self.prng = self.spawn(None)                        # replace prng



################################################################################
################################################################################
################################################################################

# private instance of RandomStateWrapper; public should not directly modify
_prng = RandomStateWrapper()

# expose variables of _prng to the public
default_BitGenerator = _prng.default_BitGenerator
seed_nbytes = _prng.seed_nbytes

# expose methods of _prng to the public
beta = _prng.prng.beta
binomial = _prng.prng.binomial
bytes = _prng.prng.bytes
chisquare = _prng.prng.chisquare
choice = _prng.prng.choice
dirichlet = _prng.prng.dirichlet
exponential = _prng.prng.exponential
f = _prng.prng.f
gamma = _prng.prng.gamma
geometric = _prng.prng.geometric
gumbel = _prng.prng.gumbel
hypergeometric = _prng.prng.hypergeometric
integers = _prng.prng.integers
laplace = _prng.prng.laplace
logistic = _prng.prng.logistic
lognormal = _prng.prng.lognormal
logseries = _prng.prng.logseries
multinomial = _prng.prng.multinomial
multivariate_hypergeometric = _prng.prng.multivariate_hypergeometric
multivariate_normal = _prng.prng.multivariate_normal
negative_binomial = _prng.prng.negative_binomial
noncentral_chisquare = _prng.prng.noncentral_chisquare
noncentral_f = _prng.prng.noncentral_f
normal = _prng.prng.normal
pareto = _prng.prng.pareto
permutation = _prng.prng.permutation
poisson = _prng.prng.poisson
power = _prng.prng.power
# rand = _prng.prng.rand
# randint = _prng.prng.randint
# randn = _prng.prng.randn
random = _prng.prng.random
# random_integers = _prng.prng.random_integers
# random_sample = _prng.prng.random_sample
rayleigh = _prng.prng.rayleigh
seed = _prng.seed
shuffle = _prng.prng.shuffle
spawn = _prng.spawn
standard_cauchy = _prng.prng.standard_cauchy
standard_exponential = _prng.prng.standard_exponential
standard_gamma = _prng.prng.standard_gamma
standard_normal = _prng.prng.standard_normal
standard_t = _prng.prng.standard_t
triangular = _prng.prng.triangular
uniform = _prng.prng.uniform
vonmises = _prng.prng.vonmises
wald = _prng.prng.wald
weibull = _prng.prng.weibull
zipf = _prng.prng.zipf
# Two legacy that are trivial wrappers around random_sample
# sample = _prng.prng.random_sample
# ranf = _prng.prng.random_sample
