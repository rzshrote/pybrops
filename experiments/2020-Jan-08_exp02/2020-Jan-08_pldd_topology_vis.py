# The purpose of this experiment is to examine the relationship between
# paa, pafd and pa as they vary with target frequency

import numpy, pandas, time, os, sys, shutil, errno
import seaborn
from matplotlib import pyplot
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

################################################################################
pldd_filename = "pldd_topology.tsv"

if not os.path.exists(pldd_filename):
    raise FileNotFoundError(
        errno.ENOENT,
        os.strerror(errno.ENOENT),
        pldd_filename
    )

pts_df = pandas.read_csv(pldd_filename, sep='\t', index_col=False)

################################################################################
# plot the figure
fig = pyplot.figure()
ax = Axes3D(fig)
ax.set_xlabel('PAA')
ax.set_ylabel('PLDD')
ax.set_zlabel('tld')

def init():
    ax.scatter(
        pts_df["paa"],
        pts_df["pldd"],
        pts_df["tld"],
        marker='.',
        s=5,
        alpha=0.5
    )
    return fig,

def animate(i):
    ax.view_init(elev=10.0, azim=i)
    print(i)
    return fig,

anim = animation.FuncAnimation(
    fig,
    animate,
    init_func=init,
    frames=360,
    interval=20,
    blit=True
)
anim.save("animation.mp4", fps=30, extra_args=['-vcodec', 'libx264'])
quit()
pyplot.clf()
