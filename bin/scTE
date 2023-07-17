#!/usr/bin/env python
import pandas as pd
import multiprocessing
from functools import partial
import logging
import os, sys, glob, datetime, time, gzip
import argparse
import collections
from math import log
sys.path.append(os.path.join(os.path.split(sys.argv[0])[0], '../'))
from scTE.miniglbase import genelist, glload, location
from scTE.annotation import annoGtf
from scTE.base import *

def prepare_parser():
    desc = "hahaha..."

    exmp = "Example: scTE <-i scRNA.sorted.bam> <-o out> [--min_genes 200] [--min_counts 400] [-p 4] <-x mm10.exclusive.idx>"

    parser = argparse.ArgumentParser(prog='scTE',description=desc, epilog=exmp)

    optional = parser._action_groups.pop()

    optional.add_argument('--min_genes', dest='genenumber',metavar='INT', type=int,default=200,
                        help='Minimum number of genes expressed required for a cell to pass filtering. Default: 200')

    optional.add_argument('--min_counts', dest='countnumber',metavar='INT', type=int,
                        help='Minimum number of counts required for a cell to pass filtering. Default: 2*min_genes')

    optional.add_argument('--expect-cells', dest='cellnumber',metavar='INT', type=int,  default=10000,
                        help='Expected number of cells. Default: 10000')

    optional.add_argument('-f','--format', metavar='input file format', dest='format', type=str, nargs='?', default='BAM', choices=['BAM','SAM'],
                        help='Input file format: BAM or SAM. DEFAULT: BAM')

    optional.add_argument('-CB', dest='CB', type=str, nargs='?', default='CR', choices=['CR','CB','False'],
                        help='Set to false to ignore for cell barcodes, it is useful for SMART-seq. If you set CB=False, it also will set UMI=False by default, Default: CR')

    optional.add_argument('-UMI', dest='UMI', type=str, nargs='?', default='UR', choices=['UR','UB','False'],
                        help='Set to false to ignore for UMI, it is useful for SMART-seq. Default: True')

    optional.add_argument('--keeptmp', dest='keeptmp', type=str, nargs='?', default='False', choices=['True','False'],
                        help='Keep the _scTEtmp file, which is useful for debugging. Default: False')

    optional.add_argument('--hdf5', dest='hdf5', type=str, nargs='?', default='False', choices=['True','False'],
                        help='Save the output as .h5ad formatted file instead of csv file. Default: False')

    optional.add_argument('-p','--thread', metavar='INT', dest='thread', type=int, default=1,
                        help='Number of threads to use, Default: 1')

    optional.add_argument('-v','--version', action='version', version='%(prog)s 1.0')

    required = parser.add_argument_group('required arguments')

    required.add_argument('-i','--input', dest='input', type=str, nargs='+', required=True,
                        help='Input file: BAM/SAM file from CellRanger or STARsolo, the file must be sorted by chromosome position')

    required.add_argument('-x', dest='annoglb',nargs='+', required=True,
                        help='The filename of the index for the reference genome annotation.')

#     required.add_argument('-g','--genome', metavar='genome', dest='genome', type=str, nargs='?', default='mm10', choices=['hg38','mm10',], required=True,
#                         help='"hg38" for human, "mm10" for mouse')

    required.add_argument('-o','--out', dest='out', nargs='?', required=True, help='Output file prefix')

    parser._action_groups.append(optional)
    optional = parser.add_argument_group('optional arguments')
    optional

    return parser

def main():
    """Start scTEs......parse options......"""

    timestart=datetime.datetime.now()
    args=read_opts(prepare_parser())

    # Fix up the UMI/CB booleans:
#     if args.UMI == 'True': args.UMI = True
#     else: args.UMI = False
#     if args.CB == 'True': args.CB = True
#     else: args.CB = False
    if args.hdf5 == 'True': args.hdf5 = True
    else: args.hdf5 = False

    info = args.info
    error = args.error

    assert sys.version_info >= (3, 6), 'Python >=3.6 is required'

    info(args.argtxt + "\n")

    outname = args.out.split('/')[-1:][0]

    info("Loading the genome annotation index... %s"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    allelement, chr_list, all_annot, glannot = Readanno(filename=outname, annoglb=args.annoglb[0]) #genome=args.genome
    print(sorted(chr_list))
    info("Finished loading the genome annotation index... %s \n"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    info("Processing BAM/SAM files ...%s"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    if len(args.input) == 1 and ',' in args.input[0]:
        args.input=args.input[0].split(',')

    if not os.path.exists('%s_scTEtmp/o1'%outname):
        os.system('mkdir -p %s_scTEtmp/o1'%outname)

    for k in args.input:
        checkCBUMI(filename=k,out=outname,CB=args.CB,UMI=args.UMI)
    info("Input SAM/BAM file appears to be valid")

    if len(args.input) > 1:
        info('Using parabam2bed as more than 1 input BAM')
        pool=multiprocessing.Pool(processes=args.thread)
        partial_work = partial(Para_bam2bed, CB=args.CB, UMI=args.UMI,out=outname)
        pool.map(partial_work, args.input)
        os.system('gunzip -c -f %s_scTEtmp/o0/*.bed.gz | gzip > %s_scTEtmp/o1/%s.bed.gz' % (outname,outname,outname))
    
    else:
        print(args.CB,args.UMI,'good\n')
        Bam2bed(args.input[0], args.CB, args.UMI, outname, args.thread)
    info("Done BAM/SAM files processing ...%s \n"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    info("Splitting ...%s"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    if args.thread == 1: #Single thread path, mainly
        # This is useful for testing optimsations, as the multiprocessing path the profile
        # Just gets locked up in {method 'acquire' of '_thread.lock' objects}
        info('Executing single thread path')
        whitelist = splitAllChrs(chr_list, filename=outname, genenumber=args.genenumber, countnumber=args.countnumber, UMI=args.UMI)
    else:
        info('Executing multiple thread path with %s threads' % args.thread)
        pool=multiprocessing.Pool(processes=args.thread)
        partial_work = partial(splitChr, filename=outname, CB=args.CB, UMI=args.UMI)
        pool.map(partial_work, chr_list)
        whitelist = filterCRs(filename=outname, genenumber=args.genenumber, countnumber=args.countnumber)

    info("Finished processing sample files %s \n"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    info("Fetching from the annotation index... %s"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    if args.thread == 1: #Single thread path
        for chrom in chr_list:
            align(chr=chrom, filename=outname, all_annot=None, glannot=glannot, whitelist=whitelist) #CB=args.CB

    else: # Multiprocessing path:
        pool = multiprocessing.Pool(processes=args.thread)
        partial_work = partial(align, filename=outname, all_annot=all_annot, glannot=None, whitelist=whitelist) # send a copy of the index,  CB=args.CB
        pool.map(partial_work, chr_list)

    if not os.path.exists('%s_scTEtmp/o4'%outname):
        os.system('mkdir -p %s_scTEtmp/o4'%outname)
    os.system('gunzip -c -f %s_scTEtmp/o3/%s.*.bed.gz | gzip > %s_scTEtmp/o4/%s.bed.gz' % (outname,outname,outname,outname))
    info("Done fetching... %s \n"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    info("Calculating expression... %s"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    len_res, genenumber, filename = Countexpression(filename=args.out, allelement=allelement, genenumber=args.genenumber, cellnumber=args.cellnumber, hdf5=args.hdf5)
    if args.hdf5 == True:
        info('Detect {0} cells expressed at least {1} genes, results output to {2}.h5ad'.format(len_res, genenumber, filename))
    else:
        info('Detect {0} cells expressed at least {1} genes, results output to {2}.csv'.format(len_res, genenumber, filename))
    
    info("Finished calculating expression %s"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    if args.keeptmp == 'True':
        pass
    else:
        os.system('rm -rf %s_scTEtmp'%outname)

    timeend = datetime.datetime.now()
    info("Done with %s\n" % timediff(timestart,timeend))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt !\n")
        sys.exit(0)


