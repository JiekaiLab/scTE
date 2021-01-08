"""
Utilities

Various utilities to support the genome scanning scripts.

MAny of these predate glbase3, but are a little tricky to remove as I am not sure where
they are used (if at all).

So excuse the terrible code in places. I will deprecate occasional functions from this.

R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
degenerate character N=[ATCG]
3-fold degenerate motifs re not used like the Lander paper.

"""

import sys, os, pickle

from . import config

def glload(filename):
    """
    **Purpose**
        Load a glbase binary file
        (Actually a Python pickle)

    **Arguments**
        filename (Required)
            the filename of the glbase binary file to load.

    **Returns**
        The glbase object previously saved as a binary file
    """
    assert os.path.exists(os.path.realpath(filename)), "File '%s' not found" % filename

    try:
        oh = open(os.path.realpath(filename), "rb")
        newl = pickle.load(oh)
        oh.close()
    except pickle.UnpicklingError:
        raise BadBinaryFileFormatError(filename)

    # Recalculate the _optimiseData for old lists, and new features
    try:
        if newl.qkeyfind:
            pass
        if "loc" in list(newl.keys()) or "tss_loc" in list(newl.keys()): # buckets are only present if a loc key is available.
            if newl.buckets: # added in 0.381, only in objects with tss_loc or loc key.
                pass
    except Exception:
        config.log.warning("Old glb format, will rebuild buckets and/or qkeyfind, consider resaving")
        newl._optimiseData()

    try:
        cons = len(newl._conditions) # expression-like object
        config.log.info("Loaded '%s' binary file with %s items, %s conditions" % (filename, len(newl), cons))
    except AttributeError:
        config.log.info("Loaded '%s' binary file with %s items" % (filename, len(newl)))
    return(newl)
