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

nindiv = 50

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

gebv = population.gebv(objcoeff = numpy.array([1.0]))
gebv_ix = gebv.argsort()
top_names = population.taxa[gebv_ix[-nindiv:]]
population.taxa_remove(top_names, invert = True)
print("selected top %s" % nindiv)

# build cross structure
cross = pybropt.popgen.Cross(
    population = population,
    varAfn = "dihybridDH",
    sparse = False,
    crossfn = "dihybridDH",
    matefn = "ctrl",
    rallocfn = "equal",
    c = 1,
    n = 10,
    s = numpy.inf,
    t = 0,
    mem = None
)
print("created cross")

# cgs = pybropt.breed.CGS(population, cross)
# print("created CGS")

# mogm = pybropt.breed.MOGM(population, cross)
# print("created MOGM")
#
# mogs = pybropt.breed.MOGS(population, cross)
# print("created MOGS")

opv = pybropt.breed.OPV(population, cross)
print("created OPV")

# pafd = pybropt.breed.PAFD(population, cross)
# print("created PAFD")
#
# pau = pybropt.breed.PAU(population, cross)
# print("created PAU")
#
# spstd = pybropt.breed.SPstd(population, cross)
# print("created SPstd")
#
# spstda = pybropt.breed.SPstdA(population, cross)
# print("created SPstdA")
#
# wgs = pybropt.breed.WGS(population, cross)
# print("created WGS")

ss = [i for i in range(nindiv)]
sspace = pybropt.algo.CategoricalSearchSpace(ss,ss,ss,ss,ss,ss,ss,ss,ss,ss)
print("created search space")

algo = pybropt.algo.StateHC(sspace)
print("created algorithm")

# cgs_sel = cgs.simulate(k = 10, objcoeff = numpy.array([1.0]), algorithm = None)
# mogm_sel = mogm.simulate(objcoeff = 1.0, algorithm = algo)
# mogs_sel = mogs.simulate(objcoeff = 1.0, algorithm = algo)
opv.simulate(
    algorithm = algo,
    bcycle = 5,
    objcoeff = numpy.array([1.0]),
    seed = 422020
)
# pafd_sel = pafd.simulate(objcoeff = numpy.array([1.0]), algorithm = algo)
# pau_sel = pau.simulate(objcoeff = numpy.array([1.0]), algorithm = algo)
# wgs_sel = wgs.simulate(k = 10, objcoeff = numpy.array([1.0]), algorithm = None)
print("ran simulations")

opv.population_history_to_csv("test_population.csv", index = None)

# print("CGS", cgs_sel)
# print("MOGM", mogm_sel)
# print("MOGS", mogs_sel)
# print("OPV", opv_sel)
# print("PAFD", pafd_sel)
# print("PAU", pau_sel)
# print("WGS", wgs_sel)
