# the purpose of this file is to host constants in a centralized location
# to avoid duplications throughout the codebase.

from .haldane import haldane2r
from .kosambi import kosambi2r

# define a lookup table for mapping functions
GMAP_MAP2R_DICT = {
    None:       haldane2r,
    'haldane':  haldane2r,
    'kosambi':  kosambi2r
}
