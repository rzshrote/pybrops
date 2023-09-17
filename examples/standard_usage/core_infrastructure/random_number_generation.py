#!/usr/bin/env python3

import random
import numpy
from pybrops.core.random.prng import seed
from pybrops.core.random.prng import spawn
from pybrops.core.random.prng import global_prng

# seed random number generators used by PyBrOpS
seed(48823)

# extract random numbers from the Python and NumPy native PRNGs, respectively
rnd1 = random.random()
rnd2 = global_prng.random()
rnd_vec = global_prng.normal(size=2)

# re-seed and compare
seed(48823)

# assert that the seeding worked
rnd1_reseed = random.random()
rnd2_reseed = global_prng.random()
rnd_vec_reseed = global_prng.normal(size=2)

# print the results
print("Seeding example:")
print(rnd1, "==", rnd1_reseed, ":", rnd1 == rnd1_reseed)
print(rnd2, "==", rnd2_reseed, ":", rnd1 == rnd1_reseed)
print(rnd_vec, "==", rnd_vec_reseed, ":", numpy.all(rnd1 == rnd1_reseed))
print()

# spawn multiple random number streams
prng_list = spawn(5)

print("PRNG spawning example:")
print(prng_list[0].normal(size=3))
print(prng_list[1].normal(size=3))
print(prng_list[2].normal(size=3))
print(prng_list[3].normal(size=3))
print(prng_list[4].normal(size=3))
