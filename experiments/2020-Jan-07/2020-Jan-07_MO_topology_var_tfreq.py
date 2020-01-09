# The purpose of this experiment is to examine the relationship between
# allele availability and allele frequency as they vary with target frequency

import numpy, pandas, time, os, sys, shutil
import seaborn
from matplotlib import pyplot
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
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
seed = 172020
n_indiv = 100
n_sel = 20
n_loci = 10000
n_pts = 10000 # number of random positions to take data on
favorable_tfreq = numpy.linspace(0.0,1.0,n_pts)
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
tfreq = None

# make empty arrays for storing scores
pts = numpy.empty((n_pts,4), dtype='float64')

for i,t in zip(range(n_pts),favorable_tfreq):
    # generate random sample
    rsel = numpy.random.choice(indices, n_sel, replace=False)

    # recalculate tfreq
    tfreq = numpy.where(coeff >= 0, t, 1.0-t)

    # score everything
    pts[i,0] = t
    pts[i,1] = objfn.paa(rsel, geno, wcoeff, tfreq, dtype=numpy.dtype("float64"))
    pts[i,2] = objfn.pad(rsel, geno, wcoeff, tfreq, dtype=numpy.dtype("float64"))
    pts[i,3] = objfn.pad_prime(rsel, geno, wcoeff, tfreq, dtype=numpy.dtype("float64"))

# make a dataframe
pts_df = pandas.DataFrame(pts, columns=["tfreq","paa","pad","pad_prime"])

# save dataframe
pts_df.to_csv("paa_pad_"+timestamp+".tsv", index=False, sep='\t')

# make figures and write them
# fig = pyplot.figure()
# ax = Axes3D(fig)
#
# def init():
#     ax.scatter(pts_df["paa"],pts_df["pad"],pts_df["tfreq"], marker='o', s=20, alpha=0.6)
#     return fig,
#
# def animate(i):
#     ax.view_init(elev=10.0, azim=i)
#     print(i)
#     return fig,
#
# anim = animation.FuncAnimation(fig, animate, init_func=init, frames=360, interval=20, blit=True)
#
# anim.save("paa_pad_tfreq_"+timestamp+"_"+"paa_vs_pad3D"+".mp4", fps=30, extra_args=['-vcodec', 'libx264'])
# pyplot.clf()

# fig = pyplot.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(pts_df["paa"],pts_df["pad"],pts_df["tfreq"])
# ax.view_init(30,30)
# pyplot.show()
# pyplot.savefig("paa_pad_tfreq_"+timestamp+"_"+"paa_vs_pad3D"+".png")
# pyplot.clf()

# paa_vs_pad_prime = seaborn.scatterplot(x="paa", y="pad_prime", data=pts_df)
# pyplot.savefig("paa_pad_tfreq_"+timestamp+"_"+"paa_vs_pad_prime3D"+".png")
# pyplot.clf()
#

fig = pyplot.figure()
ax = Axes3D(fig)
def init():
    ax.scatter(pts_df["pad"],pts_df["pad_prime"],pts_df["tfreq"], marker='+', s=10, alpha=0.6)
    return fig,
def animate(i):
    ax.view_init(elev=10.0, azim=i)
    print(i)
    return fig,
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=360, interval=20, blit=True)
anim.save("pad_pad_prime_tfreq_"+timestamp+"_"+"paa_vs_pad3D"+".mp4", fps=30, extra_args=['-vcodec', 'libx264'])
pyplot.clf()

# fig = pyplot.figure()
# ax = Axes3D(fig)
# ax.scatter(pts_df["pad"],pts_df["pad_prime"],pts_df["tfreq"])
# ax.view_init(30,30)
# pyplot.show()
# pyplot.savefig("paa_pad_prime_tfreq_"+timestamp+"_"+"paa_vs_pad3D"+".png")
# pyplot.clf()

# pad_vs_pad_prime = seaborn.scatterplot(x="pad", y="pad_prime", data=pts_df)
# pyplot.savefig("paa_pad_tfreq_"+timestamp+"_"+"pad_vs_pad_prime3D"+".png")
# pyplot.clf()
#
# # copy this script to new script
# shutil.copyfile(script_name, "paa_pad_"+timestamp+".py")
