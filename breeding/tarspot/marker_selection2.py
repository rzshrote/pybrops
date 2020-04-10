# append paths
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))

# 3rd party
import numpy
numpy.set_printoptions(threshold=sys.maxsize)

# our libraries
import pybropt.popgen

igmap = pybropt.popgen.GeneticMap.from_egmap("tarspot.McMullen_2009_US_NAM.interpolated.egmap")
print("loaded egmap")

print(igmap.chr_grp[0:10])
print(igmap.chr_start[0:10])
print(igmap.chr_stop[0:10])
print(igmap.map_pos[0:10])
print(igmap.mkr_name[0:10])
print(igmap.map_fncode[0:10])

igmap.prune(M = 0.01)

# export interpolated GeneticMap to file
igmap.to_egmap("tarspot.McMullen_2009_US_NAM.subset.egmap")
