"""

Initialise glbase, import all the libraries, set up the environment etc.

Requires:
* numpy
* matplotlib
* scipy
* sklearn
* h5py
* networkx
"""

import sys, os

#-----------------------------------------------------------------------
# Load all of the global configuration options.
try:
    from . import config
except:
    print("Error: Fatal - glbase3 is not installed correctly, cannot find my own libraries")
    print("       Is the python 'sys.path' correct?")
    sys.exit() # no raise if I can't get errors, it's surely a fatal installation problem.

# ----------------------------------------------------------------------
# Test for availability of the core non-standard libs.
# These need to be available as the subsequent load/checking is weak/non-existent.

try:
    import numpy
    config.NUMPY_AVAIL = True
except Exception:
    raise LibraryNotFoundError("Fatal - Numpy is not available or not installed")

try:
    import scipy
    config.SCIPY_AVAIL = True
except Exception:
    raise LibraryNotFoundError("Fatal - Scipy is not available or not installed")

# ----------------------------------------------------------------------
# Now import the rest of my libraries - assumes here they are available.
# If I can get config and errors then these are probably available too.

from .utils import glload
from .location import location
from .genelist import genelist

# export all of the libraries, methods and helpers.
__all__ = ["genelist",
            'config',
            "location",
            "glload",
            ]
