"""
config.py

config must be imported before any other glbase library.

"""

import logging

# -------------- Versioning data
GLBASE_VERSION = "1.1105"

# -------------- General options

SILENT = False # set this to True to silence all glbase output. Only works at startup
DEBUG = True
do_logging = True

# flags for the availability of libraries
MATPLOTLIB_AVAIL = False # required
NUMPY_AVAIL = False # required
SCIPY_AVAIL = False # required
SKLEARN_AVAIL = False # required
H5PY_AVAIL = False # Optional.
NETWORKX_AVAIL = False # optional
PYDOT_AVAIL = False # optional
NUMEXPR_AVAIL = False # Optional
PYGRAPHVIZ_AVAIL = False # Optional

# Some simple options for printing genelists
NUM_ITEMS_TO_PRINT = 3 # number of items to print by default.
PRINT_LAST_ITEM = True

# size of buckets for collide() and overlap()
# If this is changed then glload will not work correctly.
bucket_size = 10000 # in bp - tested, seems a reasonable choice.

# -------------- set up the logger here.
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-8s: %(message)s',
                    datefmt='%m-%d %H:%M'),


log = logging.getLogger('glbase3')
log.setLevel(logging.INFO)
