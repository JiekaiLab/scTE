"""

location.py

part of glbase.

This class is an internal class that implements a more convenient way to manipulate
genomic coordiantes.

TODO:
. add a 'in' code clause e.g.:
    if 1000 in location: (see if 1000 > left & < right)
    if a_location in b_location: (exectute a collide())

"""

import copy, pickle

class location:
    def __init__(self, loc=None, chr=None, left=None, right=None):
        if isinstance(loc, location):
            # It's actually already a loc.
            # I want to copy it and leave.
            self.loc = copy.copy(loc.loc)
        else:
            if loc:
                s = loc.lower().replace(",", "") # ucsc includes commas, remove them so you can cut and paste
                t = s.split(":")
                self.loc = {"chr": t[0].strip("chr").rstrip().upper(), "left":int(t[1].split("-")[0]), "right":int(t[1].split("-")[1])}
            else:
                self.loc = {"chr": str(chr).strip("chr").rstrip().upper(), "left": int(left), "right": int(right)}
        self.__update() # make sure the locstring is valid:

    def __eq__(self, other):
        if other:
            if isinstance(other, str):
                return(str(self) == str(other.replace(",", ""))) # use string comparison.

            # use a faster ? dict comparison, or throw an exception, as this item probably not a <location>
            if self.loc["chr"] == other.loc["chr"]:
                if self.loc["left"] == other.loc["left"]:
                    if self.loc["right"] == other.loc["right"]:
                        return(True)
        return(False)

    def __lt__(self, other): # deprecated in Python3
        # Make locations sortable
        if self.loc['chr'] < other.loc['chr']:
            return True
        elif self.loc['chr'] == other.loc['chr']:
            if self.loc['left'] < other.loc['left']:
                return True
            elif self.loc['left'] == other.loc['left']: # For ties
                return False
            return False
        #self.loc['chr'] > other.loc['chr']:
        return False

    def __hash__(self):
        return(hash(self._loc_string))
    
    def __deepcopy__(self, memo):
        return(pickle.loads(pickle.dumps(self, -1))) # This is 2-3x faster and presumably uses less memory
    
    def __bool__(self):
        return(True)

    def __repr__(self):
        return("<location %s>" % (self._loc_string))

    def __len__(self):
        # work out the span.
        return(max([0, self.loc["right"] - self.loc["left"]]))

    def split(self, value=None):
        # ignores the 'value' argument completely and returns a three-ple
        return( (self.loc["chr"], self.loc["left"], self.loc["right"]) )

    def __update(self):
        self._loc_string = None
        try:
            self._loc_string = "chr%s:%s-%s" % (self.loc["chr"].strip("chr"), self.loc["left"], self.loc["right"])
        except Exception: # chr possibly sets of strings ... etc.
            self._loc_string = "chr%s:%s-%s" % (self.loc["chr"], self.loc["left"], self.loc["right"])
            # I can't import my bunch of errors, as location is used in that module. So I spoof an assert
            if not self._loc_string: # failed to make a valid string...
                raise "Bad location formatting"

    def __getitem__(self, key):
        if key == "string":
            self.__update() # only update when accessed.
            return(self._loc_string)
        elif key == "dict":
            return(self.loc)
        return(self.loc[key])

    def __setitem__(self, key, value):
        self.loc[key] = value
        self.__update()

    def __str__(self):
        return(self._loc_string)

    """
    these methods below should copy the location and send a modified version back.
    """
    def expand(self, base_pairs):
        new = copy.deepcopy(self)
        new.loc["left"] -= base_pairs
        new.loc["right"] += base_pairs
        new.__update()
        return(new)

    def expandLeft(self, base_pairs):
        new = copy.deepcopy(self)
        new.loc["left"] -= base_pairs
        new.__update()
        return(new)

    def expandRight(self, base_pairs):
        new = copy.deepcopy(self)
        new.loc["right"] += base_pairs
        new.__update()
        return(new)

    def shrink(self, base_pairs):
        new = copy.deepcopy(self)
        new.loc["left"] += base_pairs
        new.loc["right"] -= base_pairs
        new.__update()
        return(new)

    def shrinkLeft(self, base_pairs):
        new = copy.deepcopy(self)
        new.loc["left"] += base_pairs
        new.__update()
        return(new)

    def shrinkRight(self, base_pairs):
        new = copy.deepcopy(self)
        new.loc["right"] -= base_pairs
        new.__update()
        return(new)

    def pointLeft(self):
        """
        get a new location at the exact left of the coordinate
        """
        new = copy.deepcopy(self)
        new.loc["right"] = new.loc["left"]
        new.__update()
        return(new)
        
    def pointRight(self):
        """
        get a new location at the exact right of the coordinate
        """
        new = copy.deepcopy(self)
        new.loc["left"] = new.loc["right"]
        new.__update()
        return(new)

    def pointify(self):
        new = copy.deepcopy(self)
        centre = (self.loc["left"] + self.loc["right"]) // 2
        new.loc = {"chr": self.loc["chr"], "left": centre, "right": centre}
        new.__update()
        return(new)

    def collide(self, loc):
        if loc["chr"] != self["chr"]:
            return(False)
        return(self.loc["right"] >= loc.loc["left"] and self.loc["left"] <= loc.loc["right"])

    def qcollide(self, loc):
        """
        **Purpose**
            perform a collision with another location object.
            This assumes you have already checked the locations are on the same chromosome.

        **Returns**
            True or False
        """
        return(self.loc["right"] >= loc.loc["left"] and self.loc["left"] <= loc.loc["right"]) # nice one-liner

    def distance(self, loc):
        """
        **Purpose**
            calculate the distance between two locations.

        **Returns**
            an integer indicating the distance, note that
            the chromosomes should be the same or it will raise an
            exception. distance() should not be used as a test for
            overlap. use collide() for that.
        """
        assert self["chr"] == loc["chr"], "chromosomes are not the same, %s vs %s" % (self, loc)
        return(self.qdistance(loc))

    def qdistance(self, loc):
        """
        (Internal)
        ignore the assert.
        """
        centreA = (self.loc["left"] + self.loc["right"]) // 2
        centreB = (loc["left"] + loc["right"]) // 2
        return(centreA - centreB)

    def __sub__(self, loc):
        """
        **Purpose**
            Allow things like:
                
            distance = locA - locB
        """
        return(self.distance(loc))

    def offset(self, base_pairs):
        """
        get a new location offset from the 5' end by n base pairs
        returns a point location.
        """
        new = copy.deepcopy(self)
        new.loc["left"] += base_pairs
        new.loc["right"] = new.loc["left"]
        new.__update()
        return(new)

    def keys(self):
        """
        Get the keys
        """
        return([i for i in self.loc])
        
if __name__ == "__main__":
    import timeit
    
    s = "a = location(loc='chr1:1000-2000').pointify()"
    t = timeit.Timer(s, "from location import location")
    print("%.2f usec/pass" % (1000000 * t.timeit(number=100000)/100000))