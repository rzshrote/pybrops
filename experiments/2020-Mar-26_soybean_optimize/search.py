# append paths
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))

# 3rd party
import numpy

# our libraries
import pybropt.algo

sspace = pybropt.algo.CategoricalSearchSpace(
    ["garbage"], ["garbage"], "Categorical Test",
    [1,2,3,4],
    [3,4,5,6]
)

print(sspace.space)
print(type(sspace.space))

print(sspace.state)
print(type(sspace.state))

print(sspace.dim_size)
print(type(sspace.dim_size))
