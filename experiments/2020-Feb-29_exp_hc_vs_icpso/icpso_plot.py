import numpy
import pandas
from sklearn.decomposition import PCA


icpso_iter_arr = numpy.loadtxt(
    "icpso_history.tsv",
    dtype = "int64",
    delimiter = '\t',
    skiprows = 1,
    usecols = 0
)
print("load iter")

icpso_score_arr = numpy.loadtxt(
    "icpso_history.tsv",
    dtype = 'float64',
    delimiter = '\t',
    skiprows = 1,
    usecols = 1
)
print("load score")

icpso_pos_arr = numpy.loadtxt(
    "icpso_history.tsv",
    dtype = 'float64',
    delimiter = '\t',
    skiprows = 1,
    usecols = [i for i in range(2,2002)]
)
print("load pos")

hc_sa_set_history = pandas.read_csv(
    "hc_sa_set_history.tsv",
    sep = "\t"
)
print("load hc_sa_set history")

print("done loading")
