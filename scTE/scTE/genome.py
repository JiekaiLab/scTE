import os,sys,gzip,time
import pybedtools
import numpy as np
from glbase3 import genelist, glload, location

def annotation(genefile,tefile,mode,genome):

    form ={'force_tsv': True, 'loc': 'location(chr=column[0], left=column[1], right=column[2])', 'annot': 3}

    if genome == 'mm':
        genefilename = 'mm10.Genes.exon'
    elif genome == 'hs':
        genefilename = 'hg38.Genes.exon'
        
    raw={}
    if '.gz' in genefile:
        o=gzip.open(genefile,'rb')
    else:
        o=open(genefile,'rU')
    for line in o:
        if '.gz' in genefile:
            line=line.decode('ascii')
        if line.startswith('#'):
            continue
        line=line.strip().split('\t')
        if line[2]=='exon' or line[2]=='UTR':
            chr = line[0]
            left = int(line[3])
            riht =  int(line[4])
            name=line[8].split('gene_name "')[1].split('";')[0]
            if name not in raw:
                raw[name] = []
            raw[name].append([chr,left,riht])
    o.close()
    
    #clean the overlapping exons
    oh=gzip.open('%s.raw.bed.gz'%genefilename,'wt')
    for k in sorted(raw):
        E=[]
        for it in raw[k]:
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
    
    clean={}
    if '.gz' in genefile:
        o=gzip.open(genefile,'rb')
    else:
        o=open(genefile,'rU')
    for line in o:
        if '.gz' in genefile:
            line=line.decode('ascii')
        if line.startswith('#'):
            continue
        if 'protein_coding' not in line and 'lincRNA' not in line:
            continue
        line=line.strip().split('\t')
        if line[2]=='exon' or line[2]=='UTR':
            chr = line[0]
            left = int(line[3])
            riht =  int(line[4])
            name=line[8].split('gene_name "')[1].split('";')[0]
            if name not in clean:
                clean[name] = []
            clean[name].append([chr,left,riht])
    o.close()
    
    oh=gzip.open('%s.clean.bed.gz'%genefilename,'wt')
    for k in sorted(clean):
        E=[]
        for it in clean[k]:
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
    
    if mode == 'exclusive':
        pybedtools.set_tempdir('./')
        gene_file = pybedtools.BedTool('%s.clean.bed.gz'%(genefilename))
        
        te_file = pybedtools.BedTool(tefile)
        
        a = te_file.intersect(gene_file,wa=True,v=True)
        oh=gzip.open('%s.exclusive.gz'%tefile,'wt')
        for n,k in enumerate(a):
            k = str(k).strip().split('\t')
            oh.write('%s\t%s\n'%(k[0],'\t'.join([ str(i) for i in k[1:]])))
        oh.close()
    
        genes = genelist('%s.raw.bed.gz'%(genefilename), format=form, gzip=True)
        TEs = genelist('%s.exclusive.gz'%tefile, format=form, gzip=True)

        all_annot = genes + TEs
        
        if genome == 'mm':
            all_annot.save('mm10.exclusive.glb')
            annot = 'mm10.exclusive.glb'
        elif genome == 'hs':
            all_annot.save('hg38.exclusive.glb')
            annot = 'hg38.exclusive.glb'
    
    elif mode == 'inclusive':
        genes = genelist('%s.raw.bed.gz'%(genefilename), format=form, gzip=True)
        if tefile.endswith('.gz'):
            TEs = genelist(tefile, format=form, gzip=True)
        else:
            TEs = genelist(tefile, format=form)
            
        
        all_annot = genes + TEs
        
        if genome == 'mm':
            all_annot.save('mm10.inclusive.glb')
            annot = 'mm10.inclusive.glb'
        elif genome == 'hs':
            all_annot.save('hg38.inclusive.glb')
            annot = 'hg38.inclusive.glb'
    return annot


# annoGtf(genefile= 'gencode.v30.annotation.gtf.gz', tefile= 'hg38.TEs.bed.gz',mode='exclusive',genome='hs')
# annoGtf(genefile= 'gencode.v30.annotation.gtf.gz', tefile= 'hg38.TEs.bed.gz',mode='inclusive',genome='hs')
# 
# annoGtf(genefile= 'gencode.vM21.annotation.gtf.gz', tefile= 'mm10.TEs.bed.gz',mode='exclusive',genome='mm')
# annoGtf(genefile= 'gencode.vM21.annotation.gtf.gz', tefile= 'mm10.TEs.bed.gz',mode='inclusive',genome='mm')





