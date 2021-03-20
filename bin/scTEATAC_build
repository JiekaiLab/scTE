#!/usr/bin/env python

import multiprocessing
from functools import partial
import logging
import os, sys, glob, datetime, time, gzip
import argparse
import collections
from math import log
sys.path.append(os.path.join(os.path.split(sys.argv[0])[0], '../'))
from scTE.miniglbase import genelist, glload, location

chr_list = [ str(k) for k in  list(range(1,50))] + ['X','Y', 'M']

def read_opts(parser):
    args = parser.parse_args()

#     if args.mode not in ['inclusive', 'exclusive'] :
#         logging.error("Counting mode %s not supported\n" % (args.mode))
#         parser.print_help()
#         sys.exit(1)
# 
#     if args.genome not in ['mm10', 'hg38'] :
#         logging.error("Counting mode %s not supported\n" % (args.genome))
#         parser.print_help()
#         sys.exit(1)

    args.info = logging.info
    return args

def genomeIndex(genome,outname):


    form={'force_tsv':True, 'loc': 'location(chr=column[0], left=column[1], right=column[2])', 'annot': 3}
    if genome.endswith('.gz'):
        genome = genelist(genome, format=form, gzip=True)
    else:
        genome = genelist(genome, format=form)
        
    genome.save('%s.idx'%outname)

def prepare_parser():
    
    desc = "Build genome annotation index for scTE"

    exmp = "Example: scTEATAC_build -g Data/TE.bed -o mm10.te"
    
    parser = argparse.ArgumentParser(prog='scTE_build',description=desc, epilog=exmp)
    
    optional = parser._action_groups.pop()
                            
    optional.add_argument('-g','--genome', metavar='genome', dest='genome',type=str, nargs='?', required=True,
                        help='Bed file of the genome window')
    
    optional.add_argument('-o','--out', dest='out', nargs='?', help='Output file prefix, Default: the genome name')
    
    required = parser.add_argument_group('required arguments')
    
    parser._action_groups.append(optional)
    optional = parser.add_argument_group('optional arguments')
    optional
    
    return parser

def main():

    timestart=datetime.datetime.now()
    args=read_opts(prepare_parser())

    assert sys.version_info >= (3, 6), 'Python >=3.6 is required'

    info = args.info

    info("Building the scTE genome annotation index... %s"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    
    genomefile=args.genome
    genomeIndex(args.genome,args.out)
    
    info("Done genome annotation index building... %s"%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt !\n")
        sys.exit(0)



