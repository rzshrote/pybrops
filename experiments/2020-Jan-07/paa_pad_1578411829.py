# The purpose of this experiment is to examine the relationship between
# allele availability and allele frequency (spoiler: they're very weakly
# correlated). This has applications for identifying an optimal point on a
# Pareto front.

import numpy, pandas, time, os, sys, shutil
import matplotlib.pyplot as pyplot
import seaborn
# hack to append into our path the parent directory for this file
sys.path.append(
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(
                os.path.realpath(__file__)
            )
        )
    )
)
# import our libraries
from pybropt import objfn

################################################################################
# define several constants
seed = 162020
n_indiv = 100
n_sel = 20
n_loci = 10000
favorable_tfreq = 0.5
n_pts = 10000 # number of random positions to take data on
indices = numpy.array([i for i in range(n_indiv)])

# get script name
script_name = os.path.basename(__file__)

# make timestamp
timestamp = str(int(time.time()))

# seed rng
numpy.random.seed(seed)

# generate binary marker data
geno = numpy.random.binomial(1, 0.1, (2,n_indiv,n_loci)).astype('uint8')

# generate marker coefficients
coeff = numpy.random.normal(0, 1, n_loci)

# calculate marker weight coefficients
wcoeff = numpy.absolute(coeff)

# make target allele frequency
tfreq = numpy.where(coeff >= 0, favorable_tfreq, 1.0-favorable_tfreq)

# make empty arrays for storing scores
pts = numpy.empty((n_pts,3), dtype='float64')

for i in range(n_pts):
    # generate random sample
    rsel = numpy.random.choice(indices, n_sel, replace=False)

    # score everything
    pts[i,0] = objfn.paa(rsel, geno, wcoeff, tfreq, dtype=numpy.dtype("float64"))
    pts[i,1] = objfn.pad(rsel, geno, wcoeff, tfreq, dtype=numpy.dtype("float64"))
    pts[i,2] = objfn.pad_prime(rsel, geno, wcoeff, tfreq, dtype=numpy.dtype("float64"))

# make a dataframe
pts_df = pandas.DataFrame(pts, columns=["paa","pad","pad_prime"])

# save dataframe
pts_df.to_csv("paa_pad_"+timestamp+".tsv", index=False, sep='\t')

# make figures and write them
paa_vs_pad = seaborn.scatterplot(x="paa", y="pad", data=pts_df)
pyplot.savefig("paa_pad_"+timestamp+"_"+"paa_vs_pad"+".png")
pyplot.clf()

paa_vs_pad_prime = seaborn.scatterplot(x="paa", y="pad_prime", data=pts_df)
pyplot.savefig("paa_pad_"+timestamp+"_"+"paa_vs_pad_prime"+".png")
pyplot.clf()

pad_vs_pad_prime = seaborn.scatterplot(x="pad", y="pad_prime", data=pts_df)
pyplot.savefig("paa_pad_"+timestamp+"_"+"pad_vs_pad_prime"+".png")
pyplot.clf()

# copy this script to new script
shutil.copyfile(script_name, "paa_pad_"+timestamp+".py")
