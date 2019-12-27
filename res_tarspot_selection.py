# import 3rd party libraries
import numpy
import time
import pandas

# import our libraries
from pybropt import gs
from pybropt import stpfn

################################################################################
# open files
tarspot_df = pandas.read_csv(
    "tarspot_selection.tsv",
    sep='\t',
    header=0,
    index_col=0
)

# grab the coefficients as numpy.ndarray
tarspot_coeff = tarspot_df.loc[:,"effect"].values.copy()

# grab the loci names as numpy.ndarray
tarspot_loci = tarspot_df.loc[:,"effect"].keys().values.copy()

# grab taxa data as numpy.ndarray
tarspot_taxa = tarspot_df.iloc[:,1:].keys().values.copy()

# convert genotypes in to (1, 323, 100) numpy.ndarray
tarspot_geno = tarspot_df.iloc[:,1:].values.T[None,:].copy()

seed = 112419
n_sel = 10

################################################################################

# calculate results
results, history = gs.opv(
    geno = tarspot_geno,
    coeff = -1*tarspot_coeff,
    sel_size = n_sel,
    algorithm = "hc_sa_set",
    algorithm_varg = None,
    #  {
    #     'n': 1000,
    #     'states': numpy.tile(numpy.arange(tarspot_geno.shape[1]), (n_sel,)),
    #     'dim_sizes': numpy.repeat(tarspot_geno.shape[1], n_sel),
    #     'inertia_wt': 0.33,
    #     'accel_coeff_pbest': 0.33,
    #     'accel_coeff_gbest': 0.34,
    #     'scale_factor': 0.1,
    #     'stpfn': lambda x: (x < 30),
    # },
    seed = None,
    nthreads = 1,
    zwidth = 3,
    verbose = True
)

print(results)
print(history)

indices = history.iloc[len(history)-1,2:].values.astype(numpy.int).copy()

print("Selection:", tarspot_taxa[indices])
