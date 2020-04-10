# append paths
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))

# 3rd party
import numpy
numpy.set_printoptions(threshold=sys.maxsize)

# our libraries
import pybropt.popgen

def load_file():
    return pybropt.popgen.GeneticMap.from_egmap(
        "tarspot.McMullen_2009_US_NAM.interpolated.egmap",
        'kosambi'
    )

# make
# igmap = load_file()
# igmap.prune(M = 0.01)
# igmap.to_egmap("tarspot.M01.egmap")
#
# igmap = load_file()
# igmap.prune(M = 0.005)
# igmap.to_egmap("tarspot.M005.egmap")
#
# igmap = load_file()
# igmap.prune(nt = 100000)
# igmap.to_egmap("tarspot.nt100000.egmap")
#
# igmap = load_file()
# igmap.prune(nt = 50000)
# igmap.to_egmap("tarspot.nt50000.egmap")
#
# igmap = load_file()
# igmap.prune(nt = 100000, M = 0.01)
# igmap.to_egmap("tarspot.nt100000.M01.egmap")
#
# igmap = load_file()
# igmap.prune(nt = 100000, M = 0.005)
# igmap.to_egmap("tarspot.nt100000.M005.egmap")
#
# igmap = load_file()
# igmap.prune(nt = 50000, M = 0.01)
# igmap.to_egmap("tarspot.nt50000.M01.egmap")
#
# igmap = load_file()
# igmap.prune(nt = 50000, M = 0.005)
# igmap.to_egmap("tarspot.nt50000.M005.egmap")

igmap = pybropt.popgen.GeneticMap.from_egmap("test.egmap", 'kosambi')
igmap.prune(nt = 100000, M = 0.01)
igmap.to_egmap("test.nt100000.M01.egmap")
