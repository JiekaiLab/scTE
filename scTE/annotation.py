import os,sys,gzip,time
import numpy as np
from scTE.miniglbase import genelist, glload, location

form ={'force_tsv': True, 'loc': 'location(chr=column[0], left=column[1], right=column[2])', 'annot': 3}

def cleanexon(filename, genefilename, exons):
    if not os.path.exists('%s_scTEtmp/index'%filename):
        os.system('mkdir -p %s_scTEtmp/index'%filename)

    oh=gzip.open('%s_scTEtmp/index/%s.bed.gz'%(filename,genefilename),'wt')
    for k in sorted(exons):
        E=[]
        for it in exons[k]:
            E+=list(range(it[1],it[2]))
        E=sorted(set(E))

        s=0
        tmp=[]
        for id in range(0,len(E)-1):
            if E[id+1]-E[id] >1:
                en=id
                tmp.append([E[s],E[en]])
                s=en+1
        tmp.append([E[s],E[id+1]])

        for item in tmp:
            oh.write('%s\t%s\t%s\t%s\n'%(it[0],item[0],item[1],k))
    oh.close()

def annoGtf(filename, genefile, tefile, mode):

    genefilename = genefile.split('/')[-1:][0].replace('.gtf','').replace('.gz','')
    tefilename = tefile.split('/')[-1:][0].replace('.bed','').replace('.gz','')

    raw = {}
    clean = {}
    if '.gz' in genefile:
        o = gzip.open(genefile,'rb')
    else:
        o=open(genefile,'rU')
    for l in o:
        if '.gz' in genefile:
            l=l.decode('ascii')
        if l.startswith('#'):
            continue
        t=l.strip().split('\t')
        if t[2]=='exon' or t[2]=='UTR':
            chr = t[0].replace('chr','')
            left = int(t[3])
            riht =  int(t[4])
            name=t[8].split('gene_name "')[1].split('";')[0]

            if name not in raw:
                raw[name] = []
            raw[name].append([chr,left,riht])

            if 'protein_coding' not in l and 'lincRNA' not in l:
                continue
            if name not in clean:
                clean[name] = []
            clean[name].append([chr,left,riht])
    o.close()

    cleanexon(filename,'%s.raw'%genefilename,raw)
    cleanexon(filename,'%s.clean'%genefilename,clean)

    if mode == 'exclusive':
        gene ={}
        o = gzip.open('%s_scTEtmp/index/%s.clean.bed.gz'%(filename,genefilename),'rb')
        for l in o:
            t = l.decode('ascii').strip().split('\t')
            chr = t[0].replace('chr','')
            left = int(t[1])
            rite = int(t[2])

            left_buck = int((left-1)/10000) * 10000
            right_buck = int((rite)/10000) * 10000
            buckets_reqd = range(left_buck, right_buck+10000, 10000)

            if chr not in gene:
                gene[chr] = {}

            if buckets_reqd:
                for buck in buckets_reqd:
                    if buck not in gene[chr]:
                        gene[chr][buck] = []
                    gene[chr][buck].append([left, rite])
        o.close()

        noverlap = []
        if '.gz' in tefile:
            o = gzip.open(tefile,'rb')
        else:
            o = open(tefile,'rU')
        for n,l in enumerate(o):
            if '.gz' in tefile:
                l = l.decode('ascii')
            t = l.strip().split('\t')
            chr = t[0]
            left = int(t[1])
            rite = int(t[2])
            
            if chr not in gene:
                noverlap.append('%s\t%s\t%s\t%s\n'%(chr,left,rite,t[3]))
                continue
            
            left_buck = int((left-1)/10000) * 10000
            right_buck = int((rite)/10000) * 10000
            buckets_reqd = range(left_buck, right_buck+10000, 10000)

            if buckets_reqd:
                i = 1
                for buck in buckets_reqd:
                    if buck not in gene[chr]:
                        pass
                    else:
                        for k in gene[chr][buck]:
                            if left < k[1] and rite > k[0]:
                                i = 0
                                break
                        if i == 0:
                            break
                if i == 1:
                    noverlap.append('%s\t%s\t%s\t%s\n'%(chr,left,rite,t[3]))

        oh = gzip.open('%s_scTEtmp/index/%s.exclusive.gz'%(filename, tefilename),'wt')
        for k in noverlap:
            oh.write(k)
        oh.close()

        genes = genelist('%s_scTEtmp/index/%s.raw.bed.gz'%(filename, genefilename), format=form, gzip=True)
        TEs = genelist('%s_scTEtmp/index/%s.exclusive.gz'%(filename, tefilename), format=form, gzip=True)
        print(genes)
        print(TEs)
        
        all_annot = genes + TEs
        all_annot.save('%s_scTEtmp/index/custome.exclusive.glb'%filename)
        annot = '%s_scTEtmp/index/custome.exclusive.glb'%filename

    elif mode == 'inclusive':
        genes = genelist('%s_scTEtmp/index/%s.raw.bed.gz'%(filename,genefilename), format=form, gzip=True)
        if tefilename.endswith('.gz'):
            TEs = genelist(tefile, format=form, gzip=True)
        else:
            TEs = genelist(tefile, format=form)

        all_annot = genes + TEs
        all_annot.save('%s_scTEtmp/index/custome.inclusive.glb'%filename)
        annot = '%s_scTEtmp/index/custome.inclusive.glb'%filename

    return annot





