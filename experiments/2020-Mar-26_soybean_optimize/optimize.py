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
import pybropt.breed
import pybropt.algo

# attempt load of large genetic map
gmap = pybropt.popgen.GeneticMap.from_egmap("Song_2016.linear.M.egmap")
print("loaded egmap")

oil_model = pybropt.model.ParametricGenomicModel.from_csv(
    "rrBLUP_models_1000.csv",
    mkr_name_ix = 0,
    coeff_ix = 1
)
print("loaded oil model")

population = pybropt.popgen.Population.from_vcf(
    fname = "Song_2016_phased_chr_1000.vcf",
    base_genetic_map = gmap,
    genomic_model = oil_model,
    kind = "linear",
    fill_value = "extrapolate",
    auto_sort = True
)
print("loaded + created population")

# build cross structure
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

opv = pybropt.breed.OPV(
    population,
    cross
)
print("created CGS")

ss = [i for i in range(100)]
sspace = pybropt.algo.CategoricalSearchSpace(
    ss,
    ss,
    ss,
    ss,
    ss,
    ss,
    ss,
    ss,
    ss,
    ss
)
print("created search space")

algo = pybropt.algo.StateHC(
    sspace
)
print("created algorithm")

sel = opv.optimize(
    objcoeff = numpy.array([1.0]),
    algorithm = algo
)
print("made selections")

algo.history_to_csv("test_algo_history.csv", index = None)

print(sel)
print(population.taxa[sel])
