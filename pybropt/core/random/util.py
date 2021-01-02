import random
from numpy.random import PCG64
from numpy.random import Generator

# default BitGenerator to use for random number generation
default_BitGenerator = PCG64

# number of bytes to use when seeding prng
prng_seed_nbytes = 8

# actual generator object to use for creating random numbers
prng = Generator(default_BitGenerator())

def spawn(n = None):
    """
    Spawn new numpy PRNG streams.

    Parameters
    ----------
    n : None, int
        Number of streams to spawn using entropy from the random module.
        If None, generate one stream. Integer must be positive or zero.

    Returns
    -------
    streams : numpy.random.Generator, list of numpy.random.Generator
        New random number generator streams.
    """
    if n is None:
        b = random.randbytes(prng_seed_nbytes)      # generate random bytes
        s = int.from_bytes(b, byteorder = "big")    # convert bytes to int
        out = Generator(default_BitGenerator(s))    # create PRNG
    elif isinstance(n, int):
        if n < 0:
            raise ValueError("'{0}' must be greater than or equal to 0: {0} = {1}".format("n", n))
        bb = [random.randbytes(prng_seed_nbytes) for _ in range(n)]
        ss = [int.from_bytes(b, byteorder = "big") for b in bb]
        out = [Generator(default_BitGenerator(s)) for s in ss]
    else:
        raise TypeError("'{0}' must be of type int".format("n"))
    return out

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
    random.seed(s)      # seed random module
    prng = spawn(None)  # replace prng
