# append paths
import sys
import os
soy_dir = os.path.dirname(os.path.realpath(__file__))   # get PyBrOpt/experiments/2020-Mar-26_soybean_optimize
exp_dir = os.path.dirname(soy_dir)                      # get PyBrOpt/experiments/
PyBrOpt_dir = os.path.dirname(exp_dir)                  # get PyBrOpt/
sys.path.append(PyBrOpt_dir)                            # append PyBrOpt

# 3rd party

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
quit()
population = pybropt.popgen.Population.from_vcf(
    fname = "Song_2016_phased_chr.vcf",
    base_genetic_map = gmap,
    kind = "linear",
    fill_value = "extrapolate"
)
print("created population")

print(population.geno.shape)
