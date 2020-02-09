import pandas

# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.realpath(__file__)
    )))
)

# import our libraries
from pybropt.popgen.gmap.gmap_class import gmap


# map = gmap.from_csv(
#     "Song_2016_gmap.csv",
#     chr_grp_ix = 4,
#     chr_start_ix = 5,
#     chr_stop_ix = 6,
#     map_pos_ix = 8,
#     map_name_ix = 1
# )
#
# map.to_gmap("Song_2016.gmap")

# map2 = gmap.from_gmap("Song_2016.gmap")
#
# map2.remove_discrepancies()
#
# map2.to_gmap("Song_2016.linear.gmap")

map3 = gmap.from_gmap("Song_2016.linear.gmap")

snp50k = pandas.read_csv("soysnp50k_pos.tsv", sep='\t', header=None)

imap = map3.interpolate(snp50k, map_name_ix = 3)

imap.to_gmap("snp50k.gmap")
