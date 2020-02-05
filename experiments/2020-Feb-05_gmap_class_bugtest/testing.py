# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.realpath(__file__)
    )))
)

# import our libraries
from pybropt.popgen.gmap.gmap_class import gmap


gmap.from_csv(
    "Song_2016_gmap.csv",
    chr_grp_ix = 4,
    chr_start_ix = 5,
    chr_stop_ix = 6,
    map_pos_ix = 8,
    map_name_ix = 1
)

gmap.to_gmap("Song_2016.gmap")
