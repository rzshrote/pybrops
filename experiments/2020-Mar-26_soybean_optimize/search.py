# append paths
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))

# 3rd party
import numpy

# our libraries
import pybropt.algo

sspace = pybropt.algo.CategoricalSearchSpace(
    None, None, "Categorical Test",
    [1,2,3,4],
    [3,4,5,6],
    [4,5,6,7,8,9,0]
)

icpso = pybropt.algo.ICPSO(
    ssize = 10,
    inertia_wt = 0.33,
    pbest_comp = 0.33,
    gbest_comp = 0.34,
    scale_factor = 0.1,
    search_space = sspace
)
