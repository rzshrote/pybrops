# append paths
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))

# 3rd party
import numpy
numpy.set_printoptions(threshold=sys.maxsize)

# our libraries
import pybropt.popgen

gmap = pybropt.popgen.GeneticMap.from_egmap("McMullen_2009_US_NAM.M.egmap")
print("loaded egmap")

population = pybropt.popgen.Population.from_vcf(
    fname = "tarspot_geno.vcf.gz",
    base_genetic_map = gmap,
    genomic_model = None,
    kind = "linear",
    fill_value = "extrapolate",
    auto_sort = True
)
print("loaded + created population")
