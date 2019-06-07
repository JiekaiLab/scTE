import os,sys,gzip,time
import pybedtools
import numpy as np
from glbase3 import genelist, glload, location

def annoGtf(genefile,tefile,mode):
    print('Start building annotation index...')
    
    genefilename = genefile.split('/')[-1:][0].replace('.gtf','')
    tefilename = tefile.split('/')[-1:][0].replace('.bed','')
    print(genefilename,tefilename)
    
    if not os.path.exists('annotation/_tmp'):
        os.system('mkdir -p annotation/_tmp')

    raw={}
    if '.gz' in genefile:
        o = gzip.open(genefile,'rb')
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
    oh=gzip.open('annotation/_tmp/%s.raw.bed.gz'%genefilename,'wt')
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
        o = gzip.open(genefile,'rb')
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
    
    oh=gzip.open('annotation/_tmp/%s.clean.bed.gz'%genefilename,'wt')
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
    
    te_file = pybedtools.BedTool(tefile)
            
    if mode == 'exclusive':
        pybedtools.set_tempdir('annotation/_tmp')
        gene_file = pybedtools.BedTool('annotation/_tmp/%s.clean.bed.gz'%genefilename)
        
        a = te_file.intersect(gene_file,wa=True,v=True)
        
        oh=gzip.open('annotation/_tmp/%s.exclusive.gz'%tefilename,'wt')
        for n,k in enumerate(a):
            k = str(k).strip().split('\t')
            oh.write('%s\t%s\n'%(k[0],'\t'.join([ str(i) for i in k[1:]])))
        oh.close()
    
        form ={'force_tsv': True, 'loc': 'location(chr=column[0], left=column[1], right=column[2])', 'annot': 3}
        genes = genelist('annotation/_tmp/%s.raw.bed.gz'%(genefilename), format=form, gzip=True)
        TEs = genelist('annotation/_tmp/%s.exclusive.gz'%tefilename, format=form, gzip=True)
        
        all_annot = genes + TEs
        all_annot.save('annotation/_tmp/custome.exclusive.glb')
        annot = 'annotation/_tmp/custome.exclusive.glb'
    
    elif mode == 'inclusive':
        genes = genelist('annotation/_tmp/%s.raw.bed.gz'%(genefilename), format=form, gzip=True)
        if tefilename.endswith('.gz'):
            TEs = genelist(tefile, format=form, gzip=True)
        else:
            TEs = genelist(tefile, format=form)

        all_annot = genes + TEs
        all_annot.save('annotation/_tmp/custome.inclusive.glb')
        annot = 'annotation/_tmp/custome.inclusive.glb'
    
    return annot

