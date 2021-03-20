"""
A special ordered list of genes.
behaves like a normal list, but each element contains a heterogenous set of data.

"""

import sys, os, csv, copy, random, pickle, re, numpy, scipy, gzip

from operator import itemgetter

from . import config
from . import utils
from .location import location
from .base_genelist import _base_genelist

class genelist(_base_genelist):
    """
    **Purpose**
        This is a class container for any arrangement of heterogenous data.

        it is good for dealing with csv/tsv files with arbitrary columns - arranging
        them by keys and then allowing cross-matching and data extraction.
        genelist is the generic implementation. Other derived classes are available for
        e.g. peaklists, genome annotations and microarray (or expression) based data.

    **Arguments**
        name (Optional)
            Will default to the name of the file if a file is loaded,
            otherwise name will be set to "Generic List" by default.
            us the name argument to give it a custom nam,e

        filename (Optional)
            the name of a file to attempt to load.

        force_tsv (Optional)
            if you specify a filename to load then
            setting this argument to True will force the file to be treated
            as a tsv (tab-separated) rather than the default csv (comma
            separated).

        loadable_list (Optional, default=None)
            If you supply a list of dicts then glbase will attempt to build a genelist out of it.

        gzip (Optional, default=False)
            The input file is gzipped

        format (Required, or use format.sniffer to let glbase guess)
            Format specifiers are a mini language that can describe any valid TSV file.

            They should be a dictionary, in the form:

            {'column_name1': 1,
            'column_name2': 2}

            where the key of the dictionary will be the key name, and the value will be the
            column number.

            location data is handled a little differently, as the values may be split across
            several columns. In this case glbase is looking for a special 'location' tag with a specific format.

            Suppose the chromosome was in column 3, the left coordinate was in column 4 and right in column 5,
            to add a genome location to our format, we would add a 'loc' key, containing this info:

            {'column_name1': 1,
            'column_name2': 2,
            "loc": "location(chr=column[3], left=column[3], right=column[5])"}

            To help deal with unusal syntax in TSV or CSV files there are a list of reserved
            key names that perform some special function:


            duplicates_key
                ?

            skiplines
                do not start loading from the file until you get to line number 'skiplines': value

            debug
                print out a debug load of the file, stopping at 'debug': X line numbers

            special
                ?

            skiptill
                do not start loading the file until you see a line starting with 'skiptill' and
                start loading the file from the next line.

            force_tsv
                forse the loader to assume the file is a TSV, rather than the defaul CSV

            gtf_decorators
                This specifies the column number that contains GTF decorators, which will be split into key:value and added to the genelist

            endwith
                Stop loading the file if you see a line that contains the value specified in endwith

            __description__
                ?

            commentlines
                Ignore lines that start with this string (e.g. 'commentlines': '#' is quite commont)

            keepifxin
                ?

            __column_must_be_used
                ?

            __ignore_empty_columns
                Ignore a column if there is no value in the column, this is for when TSVs/CSVs
                are incomplete and are missing columns on specific lines, but you don't want to have
                to sanitise the TSV/CSV, and would prefer to just fill in the blank with nothing.

            As an example, here is the full format for a complete BED file:

            {"loc": "location(chr=column[0], left=column[1], right=column[2])",
            "name": 3, "score": 4, "strand": 5, "thickStart": 6, "thickEnd": 7,
            "itemRgb": 8, "blockCount": 9, "blockSizes": 10, "blockStarts": 11,
            "force_tsv": True, "skiplines": -1, "commentlines": "#"}

            see also glbase/format.py for a list of already defined format specifiers
            that you can call using:

            gl = genelist(..., format=format.full_bed)

    """
    def __init__(self, filename=None, loadable_list=None, gzip=False, **kargs):
        # This call signature is used in a few places, so modify with care
        valig_args = ["name", "format", "force_tsv",]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError(self.__init__, k)

        self.linearData = []
        self.dataByChr = None # this is private, use get by loc.
        self.debug = False
        self.name = "Generic List"
        self.metadata = {} # container for various metadata's to extract figures from.
        self.__deathline = None # Error reporting in load_CSV()
        self.__deathindx = None

        if "format" in kargs:
            format = kargs["format"] # I expect a filename = is coming.

        if "force_tsv" in kargs and kargs["force_tsv"]:
            format["force_tsv"] = True

        if filename:
            if "format" in kargs:
                self.load(filename=filename, format=format, gzip=gzip)
            else:
                raise AssertionError('Due to excessive ambiguity the sniffing function of genelists has been removed and you now MUST provide a format argument, you can reenable this feature by specifying the sniffer: format=format.sniffer')

            config.log.info("genelist(): loaded '%s' found %s items" % (filename, len(self.linearData)))
        elif loadable_list:
            self.load_list(loadable_list)

        if "name" in kargs: # Here so it overrides anything above.
            self.name = kargs["name"]

    def load(self, filename=None, format=None, gzip=False, **kargs):
        """
        **Purpose**

        load a file into the genelist. load will attempt to load the file
        based on the filename, unless a format is specified.

        **Arguments**

        filename
            absolute filename (including path) to the actual file.
            can include path short cuts (e.g. "./", "../" etc)

        format (Optional, default = "sniffer" (ie. guess))
            format specifer, see format.py, flags.py and helpers.py and the
            documentation on how to write a valid format specifier

        **Result**

        fills the genelist with the data from the file as specified by
        the format specifier.
        """
        assert filename, "No filename specified"
        assert os.path.exists(os.path.realpath(filename)), "File %s not found" % filename

        self.path = os.path.split(os.path.realpath(filename))[0]
        self.filename = os.path.split(os.path.realpath(filename))[1]
        self.fullfilename = filename
        if self.filename.find(".") != -1:
            self.name = "".join(self.filename.split(".")[:-1])
        else:
            self.name = self.filename

        if format:
            if "special" in format: # special loads
                if format["special"] == "fasta":
                    self.linearData = utils.convertFASTAtoDict(filename=filename, gzip_input=gzip)
                    # See if I can parse names into a location?
                    try:
                        for item in self.linearData:
                            item["loc"] = location(loc=item["name"])
                    except Exception:
                        pass
                    self._optimiseData()
                    return(True)
                if format["special"] == "hmmer_tbl":
                    self.linearData = _load_hmmer_tbl(filename)
                    self._optimiseData()
                    return(True)
        else:
            raise AssertionError('Due to excessive ambiguity the sniffing function of genelists has been removed and you now MUST provide a format argument')

        csv_headers = frozenset(["csv", "xls", "tsv", "txt", "bed"])
        if filename.split(".")[-1].lower() in csv_headers: # check the last one for a csv-like header
            self.loadCSV(filename=filename, format=format, gzip=gzip, **kargs)
        elif filename.split(".")[-1] in ["glb"]:
            self = glload(filename) # will this work?
        else:
            self.loadCSV(filename=filename, format=format, gzip=gzip, **kargs)

        if "force_tsv" not in kargs and "force_tsv" not in format and len(list(self.keys())) == 1:
            config.log.warning("List contains only a single key, are you sure this is not a tsv?")

        return(True) # must have made it to one - if it fails it should trigger

    def loadCSV(self, filename=None, format=None, **kargs):
        """
        **Purpose**

        load a CSV file into the genelist

        **Arguments**

        filename
            absolute filename (including path) to the actual file.
            can include path short cuts (e.g. "./", "../" etc)

        format (Optional, default = "sniffer" (ie. guess))
            format specifer, see flags.py and helpers.py and the
            documentation on how to write a valid format specifier

        force_tsv (Optional, default=False)
            if you don't send a format specifier, but want to force the
            genelist to load as a tsv, you can set this flag to True.
            NOTE: If you send a format argument then this argument is ignored!

            As another option, you can use the sniffer_tsv format specifier

        name (Optional, Default = based on the filename)
            name of the genelist, used intermittently as labels in
            figures and the like.

        **Result**

        fills the genelist with the CSV table.
        """
        assert os.path.exists(os.path.realpath(filename)), "File %s not found" % filename

        self.path = os.path.split(os.path.realpath(filename))[0]
        self.filename = os.path.split(os.path.realpath(filename))[1]
        self.fullfilename = filename
        self.name = self.filename.split(".")[0]

        # No error wrtapping;
        self._loadCSV(filename=self.fullfilename, format=format, **kargs)

    def _loadCSV(self, **kargs):
        """
        (Internal)

        Actual loadCSV()

        This is mainly so the try/except block above doesn't get completely out of control
        and allows debug code.
        """
        assert "filename" in kargs, "No filename specified"
        assert "format" in kargs, "_loadCSV requres a format specifier"
        assert os.path.exists(kargs["filename"]), "%s file not found" % kargs["filename"]

        filename = kargs["filename"]
        format = kargs["format"]

        temp_data = []
        if 'gzip' in kargs and kargs['gzip']:
            oh = gzip.open(filename, "rt")
        else:
            oh = open(filename, "rt")

        if "force_tsv" in kargs and kargs["force_tsv"]: # force_tsv takes priority
            reader = csv.reader(oh, dialect=csv.excel_tab)
        elif "force_tsv" in format and format["force_tsv"]:
            reader = csv.reader(oh, dialect=csv.excel_tab)
        elif "dialect" in format and format["dialect"]:
            reader = csv.reader(oh, dialect=format["dialect"])
        else:
            reader = csv.reader(oh)

        if "skiplines" in format:
            skiplines = format["skiplines"]
        else:
            skiplines = 0 # skip any header row by default.

        if "skiptill" in format and format["skiptill"]:
            skiptill = format["skiptill"]
        else:
            skiptill = "Done" # This is done as truth testing fails as format["skiptill"] != None

        debug_line = 0

        for index, column in enumerate(reader):
            # This is cryptically called column, when it is actually row.\
            # there is a reason for that, it is so that in the formats it appears:
            # "refseq": column[1] # see :)
            #print index, column # debug for when all else fails!
            self.__deathline = column # For debugging purposes
            self.__deathindx = index

            if not column: # if row is completely empty, so just omit.
                continue

            # passed all the tests
            temp_data.append(self._processKey(format, column))

            #print temp_data[-1] # tells you what got loaded onto the list.
        oh.close()

        self.linearData = temp_data

        self._optimiseData()
        return(True)

    def _optimiseData(self):
        """
        (Internal)
        Call me after modifying the data to bin and build the internal tables.
        """
        self.dataByChr = None
        if not self.linearData: # list is empty, not possible to optimise anything...
            return(False)

        # Guess a loc key
        loc_key = None
        if "tss_loc" in self.linearData[0]: # always use tss_loc in preference of loc, if available
            loc_key = "tss_loc"
        elif "loc" in self.linearData[0]:
            loc_key = "loc" # Don't change this though. annotate() relies on the bucket system using tss_loc

        if "tss_loc" in self.linearData[0] and "loc" in self.linearData[0]:
            config.log.warning("List contains both 'tss_loc' and 'loc'. By default glbase will use 'tss_loc' for overlaps/collisions/annotations")

        if loc_key in self.linearData[0]: # just checking the first entry.
            self.dataByChr = {}
            self.dataByChrIndexLookBack = {}
            self.buckets = {}
            for n, item in enumerate(self.linearData): # build the chromosome quick search maps.
                chr = item[loc_key]["chr"]
                if not chr in self.dataByChr:
                    self.dataByChr[chr] = []
                    self.dataByChrIndexLookBack[chr] = []
                self.dataByChr[chr].append(item)
                self.dataByChrIndexLookBack[chr].append(n) # oh sweet, sweet dirty hack...
                # I can't remember what this look-back is for, but you
                # can use it to get the linearData index even though looking at the
                # dataByChr data It is not documented for a reason!
                # New bucket system to go in here.

                if not chr in self.buckets:
                    self.buckets[chr] = {}
                # work out the bucket(s) for the location.
                # which bucket is left and right in?
                left_buck = (item[loc_key]["left"]//config.bucket_size)*config.bucket_size
                right_buck = ((item[loc_key]["right"]+config.bucket_size)//config.bucket_size)*config.bucket_size
                buckets_reqd = list(range(left_buck, right_buck, config.bucket_size))

                #print n, item[loc_key], buckets_reqd, left_buck, right_buck, len(buckets_reqd)

                for b in buckets_reqd:
                    if not b in self.buckets[chr]:
                        self.buckets[chr][b] = []
                    self.buckets[chr][b].append(n) # use index to maintain uniqueness.

        self.qkeyfind = {}
        for index, item in enumerate(self.linearData):
            for key in item:
                if not key in self.qkeyfind:
                    self.qkeyfind[key] = {}

                try:
                    if not item[key] in self.qkeyfind[key]:
                        self.qkeyfind[key][item[key]] = []
                    self.qkeyfind[key][item[key]].append(index)
                except TypeError:
                    # The item in unhashable and cannot be added to the qkeyfind
                    # This should be pretty rare if not impossible.
                    pass
        return(True)

    def isChromosomeAvailable(self, chromosome):
        """
        you must check me before trying to access dataByChr[]
        """
        if chromosome in self.dataByChr:
            return(True)
        else:
            return(False)
        return(False)

    def _findDataByKeyLazy(self, key, value): # override????? surely find?
        """
        (Internal)

        find the first matching entry of key value pair

        This version is lazy, so I take the min() and return that item
        """
        if key in self.qkeyfind:
            if value in self.qkeyfind[key]:
                return self.linearData[min(self.qkeyfind[key][value])]
        return None # not found;

    def _findDataByKeyGreedy(self, key, value): # override????? surely finditer?
        """
        finds all - returns a list
        """
        ret = []
        item_indeces = None
        if key in self.qkeyfind:
            if value in self.qkeyfind[key]:
                item_indeces = self.qkeyfind[key][value]

        if item_indeces:
            return([self.linearData[i] for i in item_indeces])
        return(None)

    def _findAllLabelsByKey(self, key):
        """
        (Internal)
        Returns a 1D list of all the values under Key.
        Most useful for things like geneList["Symbol"]
        geneList["entrez"]
        """
        return([x[key] for x in self.linearData])

    def _findByLabel(self, key, value):
        """
        (Internal)
        key is the key to search with
        toFind is some sort value to compare

        There is a subtle problem here if the user tries to find a list or other non-hashable
        object. But then It would be kind of wierd to do that...

        (This is used in at least draw._heatmap_and_plot())

        This version is deprectaed. The official internal methods are:
        _findDataByKeyLazy|Greedy()
        """
        return(self._findDataByKeyLazy(key, value)) # not found;

    def _findByLoc(self, key, loc):
        """
        (internal)
        coords should be in formal format, chrX:int(left)-int(right)
        key is unused

        lazy.
        """
        ret = []

        if self.isChromosomeAvailable(loc["chr"]):
            for item in self.dataByChr[loc["chr"]]:
                if utils.qcollide(loc["left"], loc["right"], item[key]["left"], item[key]["right"]):
                    ret.append(item)
        return(ret)

    def saveTSV(self, filename=None, **kargs):
        """
        **Purpose**
            save the geneList or similar object as a tsv
            Note: This is not always available.
            As the geneList becomes more complex it loses the ability to be
            expressed simply as a csv-file. In those cases you must use
            the save() method to save a binary representation.

            saveTSV is *generally* 'consistent' if you can succesfully save
            as a tsv then you can reload the same list as that particular
            type. Be careful though. Behaviour like this will work fine::

                microarray.saveTSV(filename="file.csv")
                anewlist = genelist(filename="file.csv")
                anewlist # this list is now a genelist and not an expression.
                         # You must load as an expression:
                rnaseq = expression(filename="file.csv")

        **Arguments**
            filename
                filename to save, including the full path.

            key_order (Optional, default=None)
                send a list, specifying the order you'd like to write the keys
                by default saveTSV() will write to the file in an essentially
                random order. But if you want to specify the order then
                send a list of key names and it will save them in that order.

                Also, you need only specify the left side of the column.
                Any unspecified columns will be appended to the right
                hand side in a random order.

                Note that saveTSV() always saves all columns of data.

            tsv (True|False)
                save as a tsv file see also saveCSV()

            no_header (Optional, default=False)
                Don't write the first line header for this file. Usually it's the list
                of keys, then the second line will be the data. If no_header is set to False
                then the first line of the file will begin with the data.

        **Result**
            returns None.
            saves a TSV representation of the genelist.
        """
        self.saveCSV(filename, tsv=True, **kargs)

    def saveCSV(self, filename=None, no_header=False, **kargs):
        """
        **Purpose**

            save the geneList or similar object as a csv
            Note: This is not always available.
            As the geneList becomes more complex it loses the ability to be
            expressed simply as a csv-file. In those cases you must use
            the save() method to save a binary representation.

            saveCSV is guaranted to be 'consistent' if you can succesfully save
            as a csv then you can reload the same list as that particular
            type. Be careful though. Behaviour like this will work fine::

                microarray.saveCSV(filename="file.csv")
                anewlist = genelist(filename="file.csv")
                anewlist # this list is now a genelist and not a microarray.
                         # You must load as a microarry:
                amicroarry = microarray(filename="file.csv")

            saving to a csv will will blank the history, and any other meta-
            data generated about the list.

        **Arguments**

            filename
                filename to save, including the full path.

            key_order (List)
                send a list, specifying the order you'd like to write the keys
                by default saveTSV() will write to the file in an essentially
                random order. But if you want to specify the order then
                send a list of key names and it will save them in that order.

                Also, you need only specify the left side of the column.
                Any unspecified columns will be appended to the right
                hand side in a random order.

                Note that saveTSV() always saves all columns of data.

            tsv (True|False)
                save as a tsv file see also saveTSV()

            no_header (Optional, default=False)
                Don't write the first line header for this file. Usually it's the list
                of keys, then the second line will be the data. If no_header is set to False
                then the first line of the file will begin with the data.

        **Result**
            returns None.
            saves a CSV representation of the geneList.
        """
        assert filename, "No filename specified"

        oh = open(filename, "w")
        if not self.linearData: # data is empty, fail graciously.
            config.log.error("csv file '%s' is empty, no file written" % filename)
            oh.close()
            return(None)

        if "tsv" in kargs and kargs["tsv"]:
            writer = csv.writer(oh, dialect=csv.excel_tab)
        else:
            writer = csv.writer(oh)

        # work out key order and the like:
        write_keys = []
        if "key_order" in kargs:
            write_keys = kargs["key_order"]
            # now add in any missing keys to the right side of the list:
            for item in list(self.keys()):
                if item not in write_keys:
                    write_keys.append(item)
        else:
            # just selece them all:
            write_keys = list(self.keys())

        if not no_header:
            writer.writerow(write_keys) # write the header row.

        for data in self.linearData:
            line = []
            for key in write_keys: # this preserves the order of the dict.
                if key in data:
                    line.append(data[key])
                else:
                    line.append("") # a blank key, fail gracefully.
            writer.writerow(line)
        oh.close()
        config.log.info("Saved '%s'" % filename)
        return(None)

    def saveBED(self, filename=None, extra_keys=None, id=None, score=None, uniqueID=False, loc_only=False, **kargs):
        """
        **Purpose**
            save the genelist in bed format
            list must have a valid loc key

            BED format files are quite diverse in their organisation.

            The official definition is here:
            http://genome.ucsc.edu/FAQ/FAQformat

            glbase saveBED will save location in columns 0, 1 and 2.

            glbase will also enforce values in columns 3, 4 and 5. A lot of downstream
            programs require at least something in columns 3, 4 and 5. glbase will make spoof
            values for these columns. Using "+" for all strand columns and 0 in columns 4 and 5.

            You can modify the data in column 4 by sending a key name with the id= argument.
            Similarly you can modify column 5 (score) by sending a key-name to score.

            Finally if your data has a strand key that will be used for column 6.

            extra_keys will then be placed in separate columns after column 6.

        **Arguments**
            filename
                Name of the file to save the bed to

            id (Optional)
                A key in the genelist to use for the id column (column 4)

            uniqueID (Optional, default=False)
                Sometime the id column must contain a unique value. If you set this to True, glbase
                will spoof a unique id but adding <name>-n to the id column, where <name> is the
                name of the genelist and n is the n'th entry in the genelist.

                if id is set to a key value then instead of the name of the genelist the key will
                be used instead and -n appended to whatever is in that key.

            score (Optional)
                A key in the genelist to use for the score column (column 5)

            extra_keys = []
                A list of extra keys to save into the bed. These will appear after column 6 (strand).

            loc_only (Optional, default=False)
                If set to True then output a BED file containing only the loc (i.e. the frist three columns:
                'chrom, left, right')

                Note that anything in extra_keys will still be respected.

        **Returns**
            A saved bed file and None
        """
        oh = open(filename, "w")

        for index, item in enumerate(self.linearData):
            todo = ["chr%s" % str(item["loc"]["chr"]), str(item["loc"]["left"]), str(item["loc"]["right"])]
            if not loc_only:
                if uniqueID:
                    if id:
                        todo += ["%s-%s" % (str(item[id]), index)]
                    else:
                        todo += ["%s-%s" % (self.name, index)]
                else:
                    if id:
                        todo += [str(item[id])]
                    else:
                        todo += ["0"]

                if score:
                    todo += [str(item[score])]
                else:
                    todo += ["0"]

                if "strand" in item:
                    todo += [str(item["strand"])]
                else:
                    todo += ["+"]

            if extra_keys:
                todo += [str(item[k]) for k in extra_keys]

            oh.write("%s\n" % ("\t".join(todo)))

        oh.close()
        config.log.info("Saved '%s' BED file" % filename)
        return(filename)

    def sort(self, key=None, reverse=False):
        """
        Sort the data into a particular order based on key.
        This sorts the list in-place in the same style as Python.
        ie.

        ret = list.sort() - is True and not the sorted list

        list.sort() - now list contains the sorted list

        **Arguments**

        key
            the key to use to sort by.
            must be some sort of sortable item

        reverse (Optional, default=False)
            By default the list is sorted smallest to largest.
            reverse = True sorts largest to smallest.

        **Result**

        returns True if it completes.
        sorts the list IN PLACE.
        """
        assert key, "No such key '%s'" % key
        assert key in self.linearData[0], "Data does not have key '%s'" % key

        self.linearData = sorted(self.linearData, key=itemgetter(key))
        if reverse:
            self.linearData.reverse()
        self._optimiseData()
        return(True)

    def reverse(self):
        """
        reverse the order of the list, in place.

        **Arguments**

        None

        **Result**

        returns True if okay or false.
        """
        self.linearData.reverse()
        self._optimiseData() # just in case.
        return(True)

    #------------------------------ Overrides --------------------------

    def __contains__(self, item):
        """
        (Override)
        There may be some obscure failures with this item - to do with
        returned lists. IF you send a [ {} ... {} ] like object
        derived from a genelist then it will fail (but then it should?)
        but if you use slices it should be okay:
        a = genelist[0] # returns a single dict
        a = genelist[0:10] # returns a new genelist
        """
        if not self.linearData:
            return(False)

        if item in self.linearData[0]:
            return(True)
        return(False)

    def __repr__(self):
        """
        (Override)
        report the underlying representation
        """
        return("glbase.genelist")

    def __str__(self):
        """
        (Override)
        give a sensible print out.
        """
        if len(self.linearData) > config.NUM_ITEMS_TO_PRINT:
            out = []
            # welcome to perl
            for index in range(config.NUM_ITEMS_TO_PRINT):
                out.append("%s: %s" % (index, ", ".join(["%s: %s" % (key, self.linearData[index][key]) for key in self.linearData[index]])))
            out = "%s\n... truncated, showing %s/%s" % ("\n".join(out), config.NUM_ITEMS_TO_PRINT, len(self.linearData))

            if config.PRINT_LAST_ITEM:
                out = "%s\n%s" % (out, "%s: %s" % (len(self.linearData), ", ".join(["%s: %s" % (key, self.linearData[-1][key]) for key in self.linearData[-1]])))

        elif len(self.linearData) == 0:
            out = "This list is Empty"

        else: # just print first entry.
            out = []
            for index in range(len(self.linearData)):
                out.append("%s: %s" % (index, ", ".join(["%s: %s" % (key, self.linearData[index][key]) for key in self.linearData[index]])))
            out = "%s\nShowing %s/%s" % ("\n".join(out), len(self.linearData), len(self.linearData))

        return(out)

    def getColumns(self, return_keys=None):
        """
        **Purpose**
            return a new genelist only containing the columns specified in return _keys (a list)
        """
        assert isinstance(return_keys, list), "return_keys must have a list"

        newl = self.shallowcopy()
        newl.linearData = []

        for item in self.linearData:
            newl.linearData.append(dict((k, item[k]) for k in return_keys)) # hack for lack of dict comprehension
            # assumes all keys are in the dict

        newl._optimiseData()
        config.log.info("getColumns: got only the columns: %s" % ", ".join(return_keys))
        return(newl)

    def getRowsByKey(self, key=None, values=None, use_re=True, case_sensitive=True, **kargs):
        """
        **Purpose**
            extract all rows from a genelist for which the values in key are in the
            list_of_items

            You can send regular expressions and they will be
            interpreted correctly.

            NOTE that getRowsByKey() is a SLOW look up.

            If you need speed use get(), which is basically a free lookup of the list
            (even for greedy searches) but does not support regular expressions

        **Arguments**
            values (Required)
                a list of items to collect

            key (Optional, default=None)
                the key to search in.
                If None, all keys are searched.

            case_sensitive (Optional, default=True)
                Set to True to make the search case sensitive. Only works if use_re=True

            use_re (Optional, default=True)
                Unset this if you want exact matches or are getting unexpected behaviour
                from the regular expressions.

        **Returns**
            A new genelist containing only those items.
        """
        assert values, "getRowsByKey: 'values' argument cannot be None"
        if not case_sensitive:
            assert use_re, 'use_re must be True if case_sensitive is False'

        if not isinstance(values, list):
            values = [values]

        # This should be made super fast with qkeyfind

        newl = self.shallowcopy()
        newl.linearData = []

        if use_re: # re-ise the items.
            if case_sensitive:
                list_of_items = [re.compile(i) for i in values]
            else:
                list_of_items = [re.compile(i, re.IGNORECASE) for i in values]

            if not key: # split here for clarity.
                for item in self.linearData:
                    for k in item: # no key specified, check all.
                        for r in list_of_items:
                            if r.search(str(item[k])):
                                newl.linearData.append(utils.qdeepcopy(item))
            else:
                for item in self.linearData:
                    for r in list_of_items:
                        if r.search(str(item[key])): # sometimes gene names accidentally get stored as numbers
                            newl.linearData.append(utils.qdeepcopy(item))
        else: # no re.
            if not key: # split here for clarity.
                for item in self.linearData:
                    for k in item: # no key specified, check all.
                        for r in values:
                            if r == item[k]:
                                newl.linearData.append(item.copy())
            else:
                for item in self.linearData:
                    for r in values:
                        if r == item[key]:
                            newl.linearData.append(item.copy())

        if newl:
            newl._optimiseData()
        else:
            config.log.info("getRowsByKey: Found 0 items")
            return(None)

        config.log.info("getRowsByKey: Found %s items" % len(newl))
        return(newl)

    def map(self, genelist=None, peaklist=None, microarray=None, genome=None, key=None,
        greedy=True, logic="and", silent=False, **kargs):
        """
        **Purpose**
            map() merges two genelist-like objects and outputs a new genelist.

            It matches, by the key, each item that overlap in the two genelist and
            returns a new genelist only containing the matching items between the two lists.

            The new genelist will inherit from 'the right', for
            example if you have a expression-object you should perform the map in this
            order::

                result = gene_list.map(genelist=expn, key="refseq") # expn is an expresion object

            'result' will now be a expression object with all the appropriate
            methods.

            If however, you did this by mistake::

                result = expn.map(genelist=gene_list, key="refseq") # expn is an expression object

            It will still work fine, but now, trying::

                result.heatmap(filename="foo.png")

            will fail, because the result is a vanilla genelist rather than
            an expression-object as you might intend.

            Also note, this method is 'greedy' by default and and will take all matching
            entries it finds. This can be changed by setting greedy=False.

        **Arguments**
            genelist
                some sort of genelist-like object,
                examples inglude genelist, expression, genome, etc

            key
                a key in common between the two lists you can use to map
                them against each other.

            image_filename (Optional)
                save a venn diagram

            venn_proportional (Optional)
                enable experimental proportional venn diagrams.

                Note that for a proper venn, both lists should be unique for
                the key you are using to match. glbase does not check that this is
                the case. This can be useful to estimate the Venn overlap so glbase remains
                silent for non-unique lists, but can occasionally give bizarre results,
                such as negative numbers in a particular condition.

            title (Optional, default = None)
                title for the venn diagram.

            greedy (Optional, default=True)
                set to True to collect all matching entries, including duplicates. (The default
                behaviour)
                If set to False then the search finds the first matching entry only

            logic (Optional, default="and")
                a logic operator to apply to the map
                Accepted operators are:

                "and" = only keep the item if it appears in both lists
                "notright" = for each item in the right hand list, only keep it if it is NOT in the left hand list

                Be aware of the strange syntax of 'notright'. This tests the item in the right list
                and only keeps it if it is NOT in the left list.

        **Result**

            returns a new genelist-like object containing the overlapping
            objects, inheriting methods from the
            right hand side of the function.

        """
        valid_args = ("genelist", "peaklist", "microarray", "key",
            "image_filename", "title", "venn_proportional", "greedy")
        for k in kargs:
            if not k in valid_args:
                raise ArgumentError(self.map, k)

        assert logic in ("and", "notright"), "logic '%s' not supported" % logic

        if repr(genelist) == "glbase.delayedlist": # delayedlists will fail an assertion
            gene_list = genelist
        else:
            assert genome or genelist or peaklist or microarray, "map(): No valid genelist specified"

        if genelist:
            gene_list = genelist
        if peaklist:
            gene_list = peaklist
        if microarray:
            gene_list = microarray
        if genome:
            gene_list = genome

        __warning_assymetric_errs = False

        assert key, "Must specify a 'key' to map the two lists"
        #assert key in gene_list.linearData[0], "'%s' key not found in provided genelist '%s'" % (key, self.name)
        #assert key in self.linearData[0], "'%s' key not found in self '%s'" % (key, self.name)
        map_key = key

        p = progressbar(len(gene_list)) # leave as len() for compatability with delayedlists
        # speed up with a triangular search?
        newl = gene_list.shallowcopy()
        if repr(genelist) == "glbase.delayedlist": # Special exception for delayedlists, need to convert to vanilla genelist:
            newl = Genelist()
            newl.name = gene_list.name

        newl.linearData = []
        for index, item in enumerate(gene_list):
            if greedy:
                results = self._findDataByKeyGreedy(map_key, item[map_key])
            else:
                results = self._findDataByKeyLazy(map_key, item[map_key])
                if results:
                    results = [results] # coerce to a single member list to simplify code below

            if results:
                if logic == "and":
                    for r in results:
                        new_entry = utils.qdeepcopy(r) # inherit from the right
                        new_entry.update(item) # Key items inherit from the right hand side

                        # add a special case for expression objects:
                        # test that both objects are actually expression objects with _conditions and ["conditions"]:
                        # Funky syntax in case I ever derive a descendent of expression:
                        if "conditions" in item and "conditions" in r: # See if conditions in both genelists:
                            # The below line will escape the rare occasions a genelist is sent that has "conditions" but no _conditions
                            if hasattr(self, '_conditions') and hasattr(gene_list, '_conditions'): # I think we can safely assume both are bonafide expression
                                if self._conditions != gene_list._conditions: # DONT do this if the conditions are identical.
                                    new_entry["conditions"] = item["conditions"] + r["conditions"]
                                    newl._conditions = gene_list._conditions + self._conditions # will update multiple times, whoops.

                                    # only look at the err keys if I am merging the conditions
                                    if "err" in item and "err" in r:
                                        if self._conditions != gene_list._conditions: # DONT do this if the conditions are identical.
                                            new_entry["err"] = item["err"] + r["err"]
                                    elif "err" in new_entry: # Only one list has an err key, output a warning and kill it.
                                        if not __warning_assymetric_errs:
                                            __warning_assymetric_errs = True
                                            config.log.warning("map(): Only one of the two lists has an 'err' key, deleting it")
                                        del new_entry["err"]

                        newl.linearData.append(new_entry)
            elif logic == "notright":
                newl.linearData.append(utils.qdeepcopy(item)) # only inherit from the right, can't inheret from the left, as no matching map

            p.update(index)

        if "image_filename" in kargs and kargs["image_filename"]:
            if not greedy:
                config.log.warning("map(): greedy=False, this can lead to inaccurate results in the Venn diagram")

            venn_proportional = False
            if "venn_proportional" in kargs and kargs["venn_proportional"]:
                venn_proportional = True

            title = "".join(kargs["image_filename"].split(".")[:-1]) # guess a title
            if "title" in kargs:
                title = kargs["title"]

            if "." in kargs["image_filename"]:
                filename_root = "".join(kargs["image_filename"].split(".")[:-1])
            else:
                filename_root = kargs["image_filename"]

            labels = {"left": self.name, "right": gene_list.name, "title": title}
            # and modify the output and draw the venn
            self.draw._vennDiagram2(len(self)-len(newl), len(gene_list)-len(newl), len(newl),
                filename="%s_venn.png" % filename_root, proportional=venn_proportional,
                labels=labels)

        if not silent:
            if logic == "notright":
                config.log.info("map: '%s' vs '%s', using '%s' via '%s', kept: %s items" % (self.name, gene_list.name, map_key, logic, len(newl)))
            else:
                config.log.info("map: '%s' vs '%s', using '%s', found: %s items" % (self.name, gene_list.name, map_key, len(newl)))

        if len(newl.linearData):
            newl._optimiseData()
            return(newl)
        return(None)

    def pointify(self, key="loc"):
        """
        Convert all of the loc coordinates to a single point, centred
        around the middle of the coordinates

        Uses a 'loc' key as the default location to pointify

        **Arguments**

            key (default = "loc")
                specify a location key to pointify.
                defaults to loc


        **Result**

            Returns a new list with 'pointified' coordinates - the coordinates
            are now a single base pair centred around the middle of the
            coordinates.
        """
        assert key in self.linearData[0], "'%s' not in this list" % key

        newl = self.deepcopy()
        for item in newl:
            item[key] = item[key].pointify()

        newl._optimiseData()

        config.log.info("Pointified peaklist '%s'" % self.name)
        return(newl)

    def expand(self, key="loc", base_pairs=None, side="both"):
        """
        Add base_pairs to the left and right of the location specified in 'key'

        Uses a 'loc' key as the default location to pointify

        **Arguments**

            key (default = "loc")
                specify a location key to pointify.
                defaults to loc

            base_pairs (Required)
                Number of base pairs

            side (Optional, default="both")
                The side to use to expand the location by.

                "both"
                loc = chromosome, left+base_pairs, right+base_pairs

                "left"
                loc = chromosome, left+base_pairs, right

                "right"
                loc = chromosome, left, right+base_pairs

        **Result**

            Returns a new list with 'pointified' coordinates - the coordinates
            are now a single base pair centred around the middle of the
            coordinates.
        """
        assert key in self.linearData[0], "'%s' not in this list" % key

        newl = self.deepcopy()
        for item in newl:
            if side == "both":
                item[key] = item[key].expand(base_pairs)
            elif side == "left":
                item[key] = item[key].expandLeft(base_pairs)
            elif side == "right":
                item[key] = item[key].expandRight(base_pairs)

        newl._optimiseData()

        config.log.info("Expanded '%s' in genelist '%s' by %s base pairs" % (key, self.name, base_pairs))
        return(newl)

    def pointLeft(self, key="loc"):
        """
        pointify the location so that it is set to the left-most base pair.

        i.e.
        loc = (chromosome, left, left)

        **Arguments**

            key (default = "loc")
                specify a location key to pointify.
                defaults to loc

        **Result**
            Returns a new list
        """
        assert key in self.linearData[0], "'%s' not in this list" % key

        newl = self.deepcopy()
        for item in newl:
            item[key] = item[key].pointLeft()

        newl._optimiseData()

        config.log.info("pointLeft genelist %s" % (self.name))
        return(newl)

    def pointRight(self, key="loc"):
        """
        pointify the location so that it is set to the left-most base pair.

        i.e.
        loc = (chromosome, right, right)

        **Arguments**
            key (default = "loc")
                specify a location key to pointify.
                defaults to loc

        **Result**
            Returns a new list
        """
        assert key in self.linearData[0], "'%s' not in this list" % key

        newl = self.deepcopy()
        for item in newl:
            item[key] = item[key].pointRight()

        newl._optimiseData()

        config.log.info("pointRight genelist %s" % (self.name))
        return(newl)

    def _collectIdenticalKeys(self, gene_list):
        """
        (Internal)
        What it says, returns a list of valid keys in common between this list and gene_list
        """
        return(list(set(self.keys()) & set(gene_list.keys())))

    def removeDuplicatesByLoc(self, mode, key="loc", delta=200):
        """
        **Purpose**
            Remove duplicates in the list based on a location tag.
            Scans all of the loc tags and removes items that overlap.

            Performs a pointify() and expand(delta) on the locations.

            This method will keep the first encountered item of a particular duplicate.
            (i.e. it is a lazy search) As such this method is not expected to give identical
            answers:

            The results you get may be influenced by the order of the list.

            Also note that any pair of peaks can only be merged a single time.
            This stops a condition where peaks can just keep merging into larger and larger peaks
            this is particularly a problem in the 'overlap' mode.

        **Arguments**
            mode (Required)
                You must set one of two modes either:
                    'pointify_expand':
                        first pointify() (i.e. take the mid point between the left and right coords)
                        and then expand(delta) (i.e. symmetrically expand the left and right by 'delta')
                        This is the old default mode
                    'overlap':
                        use the left and right span of the coords to perform an overlap.
                        Note that the new peak will take the min(l1, l2) and the max(r1, r2) where
                        l = left, r = right and 1 = the first genelist and 2 = the second genelist.

            key (Optional, default="loc")
                the loc key to use

            delta (Optional, default=200)
                The number of base pairs to search across.

                Note that delta is NOT used when the mode = 'overlap'

        **Returns**
            A new genelist containing only unique sites
        """
        if key == 'loc':
            all_keys = list(self.keys())
            assert 'tss_loc' not in all_keys, 'removeDuplicatesByLoc will give incorrect results if your genelist contains both a loc and tss_loc key. Please delete one or the other key'

        assert mode in ('pointify_expand', 'overlap'), "mode is not one of 'pointify_expand', or 'overlap'"

        if mode == 'pointify_expand':
            mask = [0] * len(self.linearData) # keep track of masked entries
            newl = []
            for index, item in enumerate(self.linearData):
                locA = item[key].pointify().expand(delta)
                if mask[index] == 0: # only do if not already masked/searched
                    # Do a collision check

                    # work out which of the buckets required:
                    left_buck = ((locA["left"]-1-delta)//config.bucket_size)*config.bucket_size # Add an extra delta for accurate bucket spanning overlaps
                    right_buck = ((locA["right"]+1+delta)//config.bucket_size)*config.bucket_size
                    buckets_reqd = list(range(left_buck, right_buck+config.bucket_size, config.bucket_size)) # make sure to get the right spanning and left spanning sites

                    # get the ids reqd.
                    loc_ids = set()
                    if buckets_reqd:
                        for buck in buckets_reqd:
                            if locA["chr"] in self.buckets:
                                if buck in self.buckets[locA["chr"]]:
                                    loc_ids.update(self.buckets[locA["chr"]][buck]) # set = unique ids
                    # loc_ids now contains all of the indeces of the items in linearData that need checking

                    for indexB in loc_ids:
                        if indexB != index:
                            other = self.linearData[indexB]
                            locB = self.linearData[indexB][key].pointify().expand(delta)
                            if locA.qcollide(locB):
                                mask[indexB] = 1

                    mask[index] = 1
                    newl.append(item) # add in the item

            ov = self.shallowcopy()
            ov.load_list(newl)
        elif mode == 'overlap':
            # get the maximum peak size for a decent estimate of delta:
            delta = 0 # This is a pad to make sure the correct buckets are selected. It should be the maximum peak width.
            for i in self:
                delta = max([delta, len(i['loc'])])
            print('max delta = ', delta)
            if delta > 1000:
                config.log.warning("removeDuplicatesByLoc: The maximum peak size is >1000 bp, performance may be poor for removeDuplicatesByLoc()")

            mask = [0] * len(self.linearData) # keep track of masked entries
            newl = []
            for index, item in enumerate(self.linearData):
                locA = item[key]
                if mask[index] == 0: # only do if not already masked/searched
                    # Do a collision check

                    # work out which of the buckets required:
                    left_buck = ((locA["left"]-1-delta)//config.bucket_size)*config.bucket_size # Add an extra delta for accurate bucket spanning overlaps
                    right_buck = ((locA["right"]+1+delta)//config.bucket_size)*config.bucket_size
                    buckets_reqd = list(range(left_buck, right_buck+config.bucket_size, config.bucket_size)) # make sure to get the right spanning and left spanning sites

                    # get the ids reqd.
                    loc_ids = set()
                    if buckets_reqd:
                        for buck in buckets_reqd:
                            if locA["chr"] in self.buckets:
                                if buck in self.buckets[locA["chr"]]:
                                    loc_ids.update(self.buckets[locA["chr"]][buck]) # set = unique ids
                    # loc_ids now contains all of the indeces of the items in linearData that need checking

                    for indexB in loc_ids:
                        if indexB != index:
                            other = self.linearData[indexB]
                            locB = self.linearData[indexB][key]
                            if locA.qcollide(locB):
                                mask[indexB] = 1

                    mask[index] = 1
                    newl.append(item) # add in the item

            ov = self.shallowcopy()
            ov.load_list(newl)

        config.log.info("Removed %s duplicates, %s remain" % ((len(self) - len(ov)), len(ov)))
        return(ov)

    def removeDuplicates(self, key=None, **kargs):
        """
        **Purpose**
            remove the duplicates in the list and returns a new list;
            keeps the first example it finds

            This will only delete duplicates within the 'key'. For example,
            these three entries in a genelist:

            1: name: Stat3, score: 20, splicing: canonical
            2: name: Stat3, score: 30, splicing: alternate
            3: name: Nanog, score: 40, splicing: alternate

            gl = gl.removeDuplicates("name")

            will give:

            1: name: Stat3, score: 20, splicing: canonical
            3: name: Nanog, score: 40, splicing: alternate

            whilst

            gl = gl.removeDuplicates("splicing")

            will result in:

            1: name: Stat3, score: 20, splicing: canonical
            2: name: Stat3, score: 30, splicing: alternate

        **Arguments**
            key
                The key in which to make search for duplicates.

        **Returns**
            The new list with the duplicates removed.
        """
        assert key, "No key specified"
        assert key in list(self.keys()), "the key '%s' was not found in this genelist" % key

        newl = self.shallowcopy()
        newl.linearData = []
        count = 0
        p = progressbar(len(self))

        for item in self.qkeyfind[key]:
            newl.linearData.append(self.linearData[min(self.qkeyfind[key][item])]) # grab first
            # Will only apply a single item (the earliest) even if there
            # IS only one of these items.

        newl._optimiseData()

        config.log.info("removeDuplicates(): %s duplicates, list now %s items long" % (len(self) - len(newl), len(newl)))
        return(newl)

    def removeExactDuplicates(self):
        """
        **Purpose**
            removes exact duplicates where all of the keys match. Keeping the first
            found copy

        **Returns**
            The new list with the duplicates removed.
        """
        newl = self.shallowcopy()
        newl.linearData = []
        count = 0

        # TODO: Could be done with: list(map(dict, frozenset(frozenset(i.items()) for i in marked_for_deletion)))
        # Lazy at the moment, this function is very rarely used. Better to speed up *ByKey()
        unq = set()
        kord = list(self.linearData[0].keys())# fix the key order

        for item in self.linearData:
            valstr = "".join([str(item[k]) for k in kord])
            if valstr not in unq:
                unq.add(valstr)
                newl.linearData.append(item) # add first item found

        newl._optimiseData()

        config.log.info("removeExactDuplicates(): %s exact duplicates" % (len(self) - len(newl)))
        return(newl)

    def load_list(self, list_to_load, name=False):
        """
        **Purpose**
            You've generated your own [{ ... }, { ...}] like list
            (A list of dicts) and you want to either reload it into
            a genelist-like object or load it into an empty genelist.
            This is the method to do that officially.

            This method should be used with great care. Some sanity
            checking is done. But not very much.

        **Arguments**
            list_to_load
                must be a list of dicts.

            name
                Allows you to change the name of the list. By default it will keep
                the previous name.

        **Returns**
            None. This is one of the few IN PLACE methods and returns
            None.
        """
        try:
            list_to_load[0]
            i = list_to_load.__iter__()
        except TypeError:
            raise AssertionError("Type Error, the list appears not to be actually a list")

        try:
            item = list_to_load[0]
            i = [item[k] for k in item]
        except Exception:
            raise AssertionError("Type Error, the list of items appears not to contain a dictionary item")

        self.linearData = pickle.loads(pickle.dumps(list_to_load, -1)) # qdeepcopy()

        # See if we have a name:
        if name:
            self.name = name

        self._optimiseData()
        return(None)

    def find(self, value):
        """
        **Purpose**
            find value in any key in the genelist.
            Note that find is lazy and returns only the first matching item it finds.
            Use map() for more accurate mapping, or get to collect all

        **Arguments**
            value
                find this 'value' anywhere in the genelist (can be in any key). Can also be any value, number
                string etc.

        **Returns**
            either the first item it finds or False
        """
        for i in self:
            for k in i:
                if i[k] == value:
                    return(i)

        return(False)
