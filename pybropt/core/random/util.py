import random
import numpy

# default BitGenerator to use for random number generation
default_BitGenerator = numpy.random.PCG64

# number of bytes to use when seeding prng
prng_seed_nbytes = 8

# actual generator object to use for creating random numbers
prng = numpy.random.Generator(default_BitGenerator())

def seed(s = None):
    """
    Seed a new pseudo-random number generator using the default BitGenerator.
    Begins by seeding the random module, then seeds pybropt.core.random.prng,
    a numpy.random.Generator object.

    Parameters
    ----------
    s : None, int
        Seed to use for the PRNG. If None, use either system time, or random
        bits from os.urandom. This behavior depends entirely on the
        functionality of random.seed().
    """
    random.seed(s)                                          # seed random module
    b = random.randbytes(prng_seed_nbytes)                  # generate random bytes
    s = int.from_bytes(b, byteorder = "big")                # convert bytes to int
    prng = numpy.random.Generator(default_BitGenerator(s))  # replace prng
