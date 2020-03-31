# append paths
import sys
import os
soy_dir = os.path.dirname(os.path.realpath(__file__))   # get PyBrOpt/experiments/2020-Mar-26_soybean_optimize
exp_dir = os.path.dirname(soy_dir)                      # get PyBrOpt/experiments/
PyBrOpt_dir = os.path.dirname(exp_dir)                  # get PyBrOpt/
sys.path.append(PyBrOpt_dir)                            # append PyBrOpt

# 3rd party
import numpy

# our libraries
import pybropt.popgen
import pybropt.model

# attempt load of large genetic map
gmap = pybropt.popgen.GeneticMap.from_egmap("Song_2016.linear.egmap")
print("loaded egmap")

oil_model = pybropt.model.ParametricGenomicModel.from_csv(
    "rrBLUP_models.csv",
    mkr_name_ix = 0,
    coeff_ix = 1
)
print("loaded oil model")

population = pybropt.popgen.Population.from_vcf(
    fname = "Song_2016_phased_chr.vcf",
    base_genetic_map = gmap,
    genomic_model = oil_model,
    kind = "linear",
    fill_value = "extrapolate"
)
print("loaded + created population")

# sort markers in their order
population.sort()
print("sorted population")

# calculate GEBVs
gebv = population.gebv(objcoeff = numpy.array([1.0]))
# for g,i in zip(gebv, population.taxa):
#     print(g,i)

# select top 100
gebv_ix = gebv.argsort()
top100_names = population.taxa[gebv_ix[-100:]]
print(top100_names)

# remove everything except the top 100 individuals
population.taxa_remove(top100_names, invert = True)
print(population.taxa)
print(len(population.taxa))
print(population.geno.shape)
# l = 0
# for t in zip(population.marker_set.mkr_name,
#              population.marker_set.map_pos,
#              population.genomic_model.mkr_name,
#              population.genomic_model.coeff):
#     if l < 50:
#         print(*t, sep='\t')
#         l += 1
#     else:
#         break

cross = pybropt.popgen.Cross(
    population = population,
    varAfn = "2wayDH",
    sparse = False,
    crossfn = "2wayDH",
    matefn = "ctrl",
    rallocfn = "equal",
    c = 1,
    n = 1,
    s = numpy.inf,
    t = 0,
    mem = None
)
print("created cross")

print(cross.varA(numpy.array([0,1])))

print("calculated varA matrix")

print(cross.varA_2wayDH_mat)
