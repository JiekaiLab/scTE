import pandas as pd
import multiprocessing
import argparse
from functools import partial
import logging
import os, sys, glob, datetime, time, gzip
import collections
from collections import defaultdict
from math import log
from scTE.miniglbase import genelist, glload, location
from scTE.annotation import annoGtf
import subprocess

import numpy as np
import scipy
import anndata as ad

def read_opts(parser):
    args = parser.parse_args()
    if args.format == "BAM" :
        args.parser = "BAM"
    elif args.format == "SAM" :
        args.parser = "SAM"
    else :
        logging.error(f"The input file must be SAM/BAM format: {args.format} !\n" )
        sys.exit(1)

    args.error = logging.critical
    args.warn = logging.warning
    args.debug = logging.debug
    args.info = logging.info

    args.argtxt ="\n".join(("Parameter list:", \
                f"Sample = {args.out}", \
#                 "Genome = %s" % (args.genome), \
                f"Reference annotation index = {args.annoglb[0]}", \
                f"Minimum number of genes required = {args.genenumber}", \
                f"Minimum number of counts required = {args.countnumber}",\
                f"Number of threads = {args.thread}",\
    ))
    return args

# def getanno(filename, genefile, tefile, genome, mode):
#     form ={'force_tsv': True, 'loc': 'location(chr=column[0], left=column[1], right=column[2])', 'annot': 3}
# 
#     if genefile == 'default' and tefile == 'default':
#         if genome == 'mm10':
#             chr_list = ['chr'+ str(i) for i in range(1,20) ] + [ 'chrX','chrY', 'chrM' ]
#             if mode == 'exclusive':
#                 if not os.path.exists('mm10.exclusive.glb'):
#                     logging.error("Did not find the annotation index mm10.exclusive.glb, you can download it from scTE github (www....) or either give the annotation with -te and -gene option \n" )
#                     sys.exit(1)
#                 all_annot = 'mm10.exclusive.glb'
#                 allelement = set(glload(all_annot)['annot'])
# 
#             elif mode == 'inclusive':
#                 if not os.path.exists('mm10.inclusive.glb'):
#                     logging.error("Did not find the annotation index mm10.inclusive.glb, you can download it from scTE github (www....) or either give the annotation with -te and -gene option \n" )
#                     sys.exit(1)
#                 all_annot = 'mm10.inclusive.glb'
#                 allelement = set(glload(all_annot)['annot'])
# 
#         elif genome == 'hg38':
#             chr_list = ['chr'+ str(i) for i in range(1,23) ] + [ 'chrX','chrY', 'chrM' ]
#             if mode == 'exclusive':
#                 if not os.path.exists('hg38.exclusive.glb'):
#                     logging.error("Did not find the annotation index hg38.exclusive.glb, you can download it from scTE github (www....) or either give the annotation with -te and -gene option \n" )
#                     sys.exit(1)
#                 all_annot = 'hg38.exclusive.glb'
#                 allelement = set(glload(all_annot)['annot'])
# 
#             elif mode == 'inclusive':
#                 if not os.path.exists('hg38.inclusive.glb'):
#                     logging.error("Did not find the annotation index hg38.inclusive.glb, you can download it from scTE github (www....) or either give the annotation with -te and -gene option \n")
#                     sys.exit(1)
#                 all_annot = 'hg38.inclusive.glb'
#                 allelement = set(glload(all_annot)['annot'])
#     else:
#         if genome in ['hg38']:
#             chr_list = ['chr'+ str(i) for i in range(1,23) ] + [ 'chrX','chrY', 'chrM' ]
# 
#         elif genome in ['mm10']:
#             chr_list = ['chr'+ str(i) for i in range(1,20) ] + [ 'chrX','chrY', 'chrM' ]
# 
#         if not os.path.isfile(tefile) :
#             logging.error("No such file: %s !\n" %(tefile))
#             sys.exit(1)
# 
#         if not os.path.isfile(genefile) :
#             logging.error("No such file: %s !\n" % (genefile))
#             sys.exit(1)
# 
#         all_annot = annoGtf(filename, genefile=genefile, tefile=tefile, mode=mode)
#         allelement = set(glload(all_annot)['annot'])
# 
#     return(allelement,chr_list,all_annot)

def Readanno(filename, annoglb): #genome
    glannot = glload(annoglb)
    allelement = set(glannot['annot'])
#     if genome in ['mm10']:
#         chr_list = ['chr'+ str(i) for i in range(1,20) ] + [ 'chrX','chrY', 'chrM' ]
#     elif genome in ['hg38']:
#         chr_list = ['chr'+ str(i) for i in range(1,23) ] + [ 'chrX','chrY', 'chrM' ]
    
    chr_list = list(set([ k['chr'] for k in glannot['loc']])) #this is useful for costume chromsome
    return(allelement, chr_list, annoglb, glannot)

def checkCBUMI(filename,out,CB,UMI):
    if CB is not False:
        subprocess.run(f'samtools view {filename} | head -100 | grep "{CB}:Z:" | wc -l > {out}_scTEtmp/o1/testCR.txt',shell=True)
        time.sleep(2) #subprocess need take some time
        o=open(f'{out}_scTEtmp/o1/testCR.txt','rU')
        for l in o:
            l=l.strip()
            if int(l) < 100:
                logging.error(f"The input file {filename} has no cell barcodes information, plese make sure the aligner have add the cell barcode key, or set CB to False")
                sys.exit(1)

    if UMI is not False:
        subprocess.run(f'samtools view {filename} | head -100| grep "{UMI}:Z:" | wc -l > {out}_scTEtmp/o1/testUMI.txt',shell=True)
        time.sleep(2)
        o=open(f'{out}_scTEtmp/o1/testUMI.txt','rU')
        for l in o:
            l=l.strip()
            if int(l) < 100:
                logging.error(f"The input file {filename} has no {UMI}:Z information, plese make sure the aligner have add the UMI key, or set UMI to False")
                sys.exit(1)

def Bam2bed(filename, CB, UMI, out, num_threads):
    if not os.path.exists(f'{out}_scTEtmp/o1'):
        os.system(f'mkdir -p {out}_scTEtmp/o1')

    sample=filename.split('/')[-1].replace('.bam','')
    if sys.platform == 'darwin': # Mac OSX has BSD sed
        switch = '-E'
    else:
        switch = '-r'

    if UMI is False:
        if CB is False:
            # Put the sample name in the barcode slot
            os.system(f'samtools view -@ {num_threads} {filename} | awk \'{{OFS="\t"}}{{print $3,$4,$4+100,"{out}"}}\' | sed {switch} \'s/^chr//g\'| gzip -c > {out}_scTEtmp/o1/{out}.bed.gz')
        else: 
            os.system(f'samtools view -@ {num_threads} {filename} | awk \'{{OFS="\t"}}{{for(i=12;i<=NF;i++)if($i~/{CB}:Z:/)n=i}}{{print $3,$4,$4+100,$n}}\' | sed {switch} \'s/{CB}:Z://g\' | sed {switch} \'s/^chr//g\' | gzip -c > {out}_scTEtmp/o1/{out}.bed.gz')
    else:
        os.system(f'samtools view -@ {num_threads} {filename} | awk \'{{OFS="\t"}}{{for(i=12;i<=NF;i++)if($i~/{CB}:Z:/)n=i}}{{for(i=12;i<=NF;i++)if($i~/{UMI}:Z:/)m=i}}{{print $3,$4,$4+100,$n,$m}}\' | sed {switch} \'s/{CB}:Z://g\' | sed {switch} \'s/{UMI}:Z://g\'| sed {switch} \'s/^chr//g\' | awk \'!x[$4$5]++\' | gzip -c > {out}_scTEtmp/o1/{out}.bed.gz')


def Para_bam2bed(filename, CB, UMI, out):
    if not os.path.exists(f'{out}_scTEtmp/o0'):
        os.system(f'mkdir -p {out}_scTEtmp/o0')

    sample=filename.split('/')[-1].replace('.bam','')
    
    if sys.platform == 'darwin': # Mac OSX has BSD sed
        switch = '-E'
    else:
        switch = '-r'
    
    print(UMI,CB)
    if UMI is False:
        if CB is False:
            os.system(f'samtools view {filename} | awk \'{{OFS="\t"}}{{print $3,$4,$4+100,"{sample}"}}\' | sed {switch} \'s/^chr//g\' | gzip > {out}_scTEtmp/o0/{sample}.bed.gz')
        else:
            os.system(f'samtools view {filename} | awk \'{{OFS="\t"}}{{for(i=12;i<=NF;i++)if($i~/{CB}:Z:/)n=i}}{{print $3,$4,$4+100,$n,$m}}\' | sed {switch} \'s/{CB}:Z://g\' | sed {switch} \'s/^chr//g\' | gzip > {out}_scTEtmp/o0/{sample}.bed.gz')
    else:
        os.system(f'samtools view {filename} | awk \'{{OFS="\t"}}{{for(i=12;i<=NF;i++)if($i~/{CB}:Z:/)n=i}}{{for(i=12;i<=NF;i++)if($i~/{UMI}:Z:/)m=i}}{{print $3,$4,$4+100,$n,$m}}\' | sed {switch} \'s/{CB}:Z://g\' | sed {switch} \'s/{UMI}:Z://g\' | sed {switch} \'s/^chr//g\' | awk \'!x[$4$5]++\' | gzip > {out}_scTEtmp/o0/{sample}.bed.gz')


def splitAllChrs(chromosome_list, filename, genenumber, countnumber, UMI=True):
    '''
    **Purpose**
        Split the data into separate beds, and count up all the times each barcode appears

        This variant uses more memory, but does it all at the same time and gets the filtered whitelist for free

    **Arguments**
        chromosome_list
            List of chromosome names

        filename (Required)
            filename stub to use for tmp files

        genenumber (Required)
            Minimum number of genes expressed required for a cell to pass filtering

        countnumber (Required)
            Minimum number of counts required for a cell to pass filtering.

        UMI (optional, default=True)
            use the UMI

    **Returns**
        The barcode whitelist
    '''

    if not os.path.exists(f'{filename}_scTEtmp/o2'):
        os.system(f'mkdir -p {filename}_scTEtmp/o2')

    chromosome_list = set([c.replace('chr', '') for c in chromosome_list])

    file_handle_in = gzip.open(f'{filename}_scTEtmp/o1/{filename}.bed.gz', 'rt')
    file_handles_out = {chr: gzip.open(f'{filename}_scTEtmp/o2/{filename}.chr{chr}.bed.gz', 'wt') for chr in chromosome_list}

    CRs = defaultdict(int)

    if UMI:
        uniques = {chrom: set([]) for chrom in chromosome_list}

    # Make a BED for each chromosome
    for line in file_handle_in:
        t = line.strip().split('\t')
        chrom = t[0].replace('chr', '') # strip chr

        if chrom not in chromosome_list: # remove the unusual chromosomes
            # Force chrMT -> chrM
            if chrom == 'MT':
                chrom = 'M'
            else:
                continue

        if UMI:
            if line in uniques[chrom]:
                continue
            uniques[chrom].add(line)
            CRs[t[3]] += 1
        else:
            CRs[t[3]] += 1

        file_handles_out[chrom].write(line)

    [file_handles_out[k].close() for k in file_handles_out]
    file_handle_in.close()

    # Because this does it all in one go, you can just filter the whitelist here now, and don't need the .count. file;
    sortcb = sorted(CRs.items(), key=lambda item:item[1], reverse=True) # Sorts by the count;

    if not countnumber:
        mincounts = 2 * genenumber
    else:
        mincounts = countnumber

    whitelist = []
    for n, k in enumerate(sortcb):
        if k[1] < mincounts:
            break
        whitelist.append(k[0])

    return set(whitelist)

def filterCRs(filename, genenumber, countnumber):
    CRs = defaultdict(int)
    for f in glob.glob(f'{filename}_scTEtmp/o2/{filename}*.count.gz'):
        o = gzip.open(f,'rt')
        for l in o:
            t = l.strip().split('\t')
            CRs[t[0]] += int(t[1])
        o.close()

    sortcb=sorted(CRs.items(),key=lambda item:item[1],reverse=True)

    if not countnumber:
        mincounts = 2* genenumber
    else:
        mincounts = countnumber

    whitelist=[]
    for n,k in enumerate(sortcb):
        if k[1] < mincounts:
            break
        whitelist.append(k[0])

    #print(len(whitelist))
    return set(whitelist)

def splitChr(chr, filename):
    if not os.path.exists(f'{filename}_scTEtmp/o2'):
        os.system(f'mkdir -p {filename}_scTEtmp/o2')

    chr=chr.replace('chr','')
    
    if chr in ['1','2','3']:
        os.system(f'gunzip -c -f {filename}_scTEtmp/o1/{filename}.bed.gz | grep -v ^{chr}\'[0-9]\' | grep ^{chr} | gzip -c > {filename}_scTEtmp/o2/{filename}.chr{chr}.bed.gz')
    else:
        os.system(f'gunzip -c -f {filename}_scTEtmp/o1/{filename}.bed.gz | grep ^{chr} | gzip -c > {filename}_scTEtmp/o2/{filename}.chr{chr}.bed.gz')

    CRs = defaultdict(int)
    o = gzip.open(f'{filename}_scTEtmp/o2/{filename}.chr{chr}.bed.gz','rt')
    for l in o:
        t = l.strip().split('\t')
        CRs[t[3]] += 1
    o.close()

    o = gzip.open(f'{filename}_scTEtmp/o2/{filename}.chr{chr}.count.gz','wt')
    for k in CRs:
        o.write(f'{k}\t{CRs[k]}\n')
    o.close()

def align(chr, filename, all_annot, glannot, whitelist): #CB
    '''
    **Purpose**
        For each read, align it to the index and assign a TE, gene.

    This is the speed critical part.

    '''
    s1 = time.time()
    chr = 'chr' + chr

    if not os.path.exists(f'{filename}_scTEtmp/o3'):
        os.system(f'mkdir -p {filename}_scTEtmp/o3')

    if not glannot: # Load separately for the multicore pipeline, share the index for the single core pipeline
        glannot = glload(all_annot)

    # Only keep the glbase parts we need.
    buckets = glannot.buckets[chr.replace('chr', '')]
    all_annot = glannot.linearData

    oh = gzip.open(f'{filename}_scTEtmp/o2/{filename}.{chr}.bed.gz', 'rt')
    res = {}
    for line in oh:
        t = line.strip().split('\t')
        barcode = t[3]
        if barcode not in whitelist:
            continue
        if barcode not in res:
            res[barcode] = defaultdict(int)

        #chrom = t[0].replace('chr', '') # Don't need as each align is already split for each chrom;
        left = int(t[1])
        rite = int(t[2])

        #loc = location(chr=chrom, left=left, right=rite)
        left_buck = ((left-1)//10000) * 10000
        right_buck = ((rite)//10000) * 10000
        buckets_reqd = range(left_buck, right_buck+10000, 10000)

        if buckets_reqd:
            loc_ids = set()
            loc_ids_update = loc_ids.update

            # get the ids reqd.
            [loc_ids_update(buckets[buck]) for buck in buckets_reqd if buck in buckets]

            result = [all_annot[index]['annot'] for index in loc_ids if (rite >= all_annot[index]['loc'].loc['left'] and left <= all_annot[index]['loc'].loc["right"])]

            if result:
                for gene in result:
                    res[barcode][gene] += 1

    oh.close()

    oh = gzip.open(f'{filename}_scTEtmp/o3/{filename}.{chr}.bed.gz', 'wt')
    for bc in sorted(res):
        for gene in sorted(res[bc]):
            oh.write('%s\t%s\t%s\n' % (bc, gene, res[bc][gene]))
    oh.close()

def Countexpression(filename, allelement, genenumber, cellnumber, hdf5):
    gene_seen = allelement

    whitelist={}
    o = gzip.open(f'{filename}_scTEtmp/o4/{filename}.bed.gz', 'rt')
    for n,l in enumerate(o):
        t = l.strip().split('\t')
        if t[0] not in whitelist:
            whitelist[t[0]] = 0
        whitelist[t[0]] += 1
    o.close()

    CRlist = []
    sortcb = sorted(whitelist.items(), key=lambda item:item[1], reverse=True)
    for n,k in enumerate(sortcb):
        if k[1] < genenumber:
            break
        if n >= cellnumber:
            break
        CRlist.append(k[0])
    CRlist = set(CRlist)

    res = {}
    genes_oh = gzip.open(f'{filename}_scTEtmp/o4/{filename}.bed.gz', 'rt')
    for n, l in enumerate(genes_oh):
        t = l.strip().split('\t')
        if t[0] not in CRlist:
            continue
        if t[0] not in res:
            res[t[0]] = {}
        if t[1] not in res[t[0]]:
            res[t[0]][t[1]] = 0
        res[t[0]][t[1]] += int(t[2])

    genes_oh.close()

    s=time.time()

    # Save out the final file

    gene_seen = list(gene_seen) # Do the sort once;
    gene_seen.sort()

    #==== save results =====
    if not hdf5: # save as csv
        res_oh = open(f'{filename}.csv', 'w')
        res_oh.write('barcodes,')
        res_oh.write('%s\n' % (','.join([str(i) for i in gene_seen])))

        for k in sorted(res):
            l = ["0"] * len(gene_seen) # Avoid all the appends
            for idx, gene in enumerate(gene_seen):
                if gene in res[k]:
                    l[idx] = str(res[k][gene])
            res_oh.write('%s,%s\n' % (k, ','.join(l)))
        res_oh.close()
    
    else: # save as hdf5
        data = []
        CBs = []
        for k in sorted(res):
            l = ["0"] * len(gene_seen) # Avoid all the appends
            for idx, gene in enumerate(gene_seen):
                if gene in res[k]:
                    l[idx] = str(res[k][gene])
            data.append(l)
            CBs.append(k)

        obs = pd.DataFrame(index = CBs)
        var = pd.DataFrame(index = gene_seen)
        adata = ad.AnnData(np.asarray(data),var = var,obs = obs)
        adata.X = scipy.sparse.csr_matrix(adata.X)
        adata.write(f'{filename}.h5ad')
    
    #========================


    return len(res), genenumber, filename

def timediff(timestart, timestop):
        t  = (timestop-timestart)
        time_day = t.days
        s_time = t.seconds
        ms_time = t.microseconds / 1000000
        usedtime = int(s_time + ms_time)
        time_hour = int(usedtime / 60 / 60 )
        time_minute = int((usedtime - time_hour * 3600 ) / 60 )
        time_second =  int(usedtime - time_hour * 3600 - time_minute * 60 )
        retstr = "%dd %dh %dm %ds"  %(time_day, time_hour, time_minute, time_second,)
        return retstr
