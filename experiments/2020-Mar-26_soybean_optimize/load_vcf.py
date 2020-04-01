# append paths
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))

# 3rd party
import numpy
numpy.set_printoptions(threshold=sys.maxsize)

# our libraries
import pybropt.popgen
import pybropt.model

# attempt load of large genetic map
gmap = pybropt.popgen.GeneticMap.from_egmap("Song_2016.linear.M.egmap")
print("loaded egmap")

oil_model = pybropt.model.ParametricGenomicModel.from_csv(
    "rrBLUP_models_1000.csv",
    mkr_name_ix = 0,
    coeff_ix = [1,2,3]
)
print("loaded oil model")

population = pybropt.popgen.Population.from_vcf(
    fname = "Song_2016_phased_chr_1000.vcf",
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
gebv = population.gebv(objcoeff = numpy.array([0.333,0.333,0.334]))
# for g,i in zip(gebv, population.taxa):
#     print(g,i)

# select top 10
gebv_ix = gebv.argsort()
top10_names = population.taxa[gebv_ix[-10:]]

# remove everything except the top 10 individuals
population.taxa_remove(top10_names, invert = True)
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
    varAfn = "dihybridDH",
    sparse = False,
    crossfn = "dihybridDH",
    matefn = "ctrl",
    rallocfn = "equal",
    c = 1,
    n = 5,
    s = numpy.inf,
    t = 0,
    mem = None
)
print("created cross")

# print(population.marker_set.map_pos)

# population.geno[:,[1,3],:] += 4

new_population = cross.mate(
    numpy.array([1,2,3,4])
)

print(len(new_population))
# print(new_population.geno)
