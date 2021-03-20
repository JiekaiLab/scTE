"""

Pack the scRNA-seq data using scanpy, prep for scran normalisation

"""

import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from rpy2.robjects.packages import importr
plt.rcParams['figure.figsize'] = (8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 10
sc.settings.autoshow = False

def sparsify(filename):
    data = pd.read_csv(filename, index_col=0, header=0)
    genes = data.columns
    cells = data.index
    data = sp.sparse.csr_matrix(data.to_numpy())
    data.astype('float32')

    '''
    oh = open('gene_names.{0}.tsv'.format(os.path.split(filename)[1]), 'w')
    for g in genes:
        oh.write('%s\n' % g)
    oh.close()
    '''

    print('Loaded {0}'.format(filename))
    ad = AnnData(data, obs={'obs_names': cells}, var={'var_names': genes})
    del data
    return ad

sam1 = sparsify("../scte_data/ss.gastrulation_E6.5_Sam1.csv.gz")    ; sam1.obs['stage'] = "E6.5"   ; sam1.obs['replicate'] = "E6.5-1"
sam2 = sparsify("../scte_data/ss.gastrulation_E6.5_Sam5.csv.gz")    ; sam2.obs['stage'] = "E6.5"   ; sam2.obs['replicate'] = "E6.5-2"
#sam3 = sparsify("../scte_data/ss.gastrulation_E6.5_Sam18.csv.gz")   ; sam3.obs['stage'] = "E6.5"   ; sam3.obs['replicate'] = "E6.5-3"
#sam4 = sparsify("../scte_data/ss.gastrulation_E6.75_Sam7.csv.gz")   ; sam4.obs['stage'] = "E6.75"  ; sam4.obs['replicate'] = "E6.75-1"
sam5 = sparsify("../scte_data/ss.gastrulation_E7.0_Sam10.csv.gz")   ; sam5.obs['stage'] = "E7.0"   ; sam5.obs['replicate'] = "E7.0-1"
#sam6 = sparsify("../scte_data/ss.gastrulation_E7.0_Sam15.csv.gz")   ; sam6.obs['stage'] = "E7.0"   ; sam6.obs['replicate'] = "E7.0-3"
sam7 = sparsify("../scte_data/ss.gastrulation_E7.0_Sam30.csv.gz")   ; sam7.obs['stage'] = "E7.0"   ; sam7.obs['replicate'] = "E7.0-4"
sam8 = sparsify("../scte_data/ss.gastrulation_E7.0_Sam31.csv.gz")   ; sam8.obs['stage'] = "E7.0"   ; sam8.obs['replicate'] = "E7.0-5"
sam9 = sparsify("../scte_data/ss.gastrulation_E7.0_Sam32.csv.gz")   ; sam9.obs['stage'] = "E7.0"   ; sam9.obs['replicate'] = "E7.0-6"
sam10 = sparsify("../scte_data/ss.gastrulation_E7.25_Sam23.csv.gz") ; sam10.obs['stage'] = "E7.25" ; sam10.obs['replicate'] = "E7.25-2"
sam11 = sparsify("../scte_data/ss.gastrulation_E7.25_Sam26.csv.gz") ; sam11.obs['stage'] = "E7.25" ; sam11.obs['replicate'] = "E7.25-3"
sam12 = sparsify("../scte_data/ss.gastrulation_E7.25_Sam27.csv.gz") ; sam12.obs['stage'] = "E7.25" ; sam12.obs['replicate'] = "E7.25-4"
sam13 = sparsify("../scte_data/ss.gastrulation_E7.5_Sam2.csv.gz")   ; sam13.obs['stage'] = "E7.5"  ; sam13.obs['replicate'] = "E7.5-1"
sam14 = sparsify("../scte_data/ss.gastrulation_E7.5_Sam3.csv.gz")   ; sam14.obs['stage'] = "E7.5"  ; sam14.obs['replicate'] = "E7.5-2"
sam15 = sparsify("../scte_data/ss.gastrulation_E7.5_Sam4.csv.gz")   ; sam15.obs['stage'] = "E7.5"  ; sam15.obs['replicate'] = "E7.5-3"
sam16 = sparsify("../scte_data/ss.gastrulation_E7.5_Sam6.csv.gz")   ; sam16.obs['stage'] = "E7.5"  ; sam16.obs['replicate'] = "E7.5-4"
sam17 = sparsify("../scte_data/ss.gastrulation_E7.5_Sam19.csv.gz")  ; sam17.obs['stage'] = "E7.5"  ; sam17.obs['replicate'] = "E7.5-5"
sam18 = sparsify("../scte_data/ss.gastrulation_E7.5_Sam20.csv.gz")  ; sam18.obs['stage'] = "E7.5"  ; sam18.obs['replicate'] = "E7.5-6"
sam19 = sparsify("../scte_data/ss.gastrulation_E7.75_Sam8.csv.gz")  ; sam19.obs['stage'] = "E7.75" ; sam19.obs['replicate'] = "E7.75-1"
sam20 = sparsify("../scte_data/ss.gastrulation_E7.75_Sam9.csv.gz")  ; sam20.obs['stage'] = "E7.75" ; sam20.obs['replicate'] = "E7.75-2"
sam21 = sparsify("../scte_data/ss.gastrulation_E7.75_Sam12.csv.gz") ; sam21.obs['stage'] = "E7.75" ; sam21.obs['replicate'] = "E7.75-3"
sam22 = sparsify("../scte_data/ss.gastrulation_E7.75_Sam13.csv.gz") ; sam22.obs['stage'] = "E7.75" ; sam22.obs['replicate'] = "E7.75-4"
sam23 = sparsify("../scte_data/ss.gastrulation_E8.0_Sam16.csv.gz")  ; sam23.obs['stage'] = "E8.0"  ; sam23.obs['replicate'] = "E8.0-1"
sam24 = sparsify("../scte_data/ss.gastrulation_E8.0_Sam33.csv.gz")  ; sam24.obs['stage'] = "E8.0"  ; sam24.obs['replicate'] = "E8.0-2"
sam25 = sparsify("../scte_data/ss.gastrulation_E8.0_Sam34.csv.gz")  ; sam25.obs['stage'] = "E8.0"  ; sam25.obs['replicate'] = "E8.0-3"
sam26 = sparsify("../scte_data/ss.gastrulation_E8.0_Sam35.csv.gz")  ; sam26.obs['stage'] = "E8.0"  ; sam26.obs['replicate'] = "E8.0-4"
sam27 = sparsify("../scte_data/ss.gastrulation_E8.25_Sam24.csv.gz") ; sam27.obs['stage'] = "E8.25" ; sam27.obs['replicate'] = "E8.25-1"
sam28 = sparsify("../scte_data/ss.gastrulation_E8.25_Sam25.csv.gz") ; sam28.obs['stage'] = "E8.25" ; sam28.obs['replicate'] = "E8.25-2"
sam29 = sparsify("../scte_data/ss.gastrulation_E8.25_Sam28.csv.gz") ; sam29.obs['stage'] = "E8.25" ; sam29.obs['replicate'] = "E8.25-3"
sam30 = sparsify("../scte_data/ss.gastrulation_E8.5_Sam17.csv.gz")  ; sam30.obs['stage'] = "E8.5"  ; sam30.obs['replicate'] = "E8.5-1"
sam31 = sparsify("../scte_data/ss.gastrulation_E8.5_Sam29.csv.gz")  ; sam31.obs['stage'] = "E8.5"  ; sam31.obs['replicate'] = "E8.5-2"
sam32 = sparsify("../scte_data/ss.gastrulation_E8.5_Sam36.csv.gz")  ; sam32.obs['stage'] = "E8.5"  ; sam32.obs['replicate'] = "E8.5-3"
sam33 = sparsify("../scte_data/ss.gastrulation_E8.5_Sam37.csv.gz")  ; sam33.obs['stage'] = "E8.5"  ; sam33.obs['replicate'] = "E8.5-4"
sam34 = sparsify("../scte_data/ss.gastrulation_mixed_Sam21.csv.gz") ; sam34.obs['stage'] = "mixed" ; sam34.obs['replicate'] = "mixed-1"
sam35 = sparsify("../scte_data/ss.gastrulation_mixed_Sam22.csv.gz") ; sam35.obs['stage'] = "mixed" ; sam35.obs['replicate'] = "mixed-2"

print('Loaded Samples...')

# Do very simple prefiltering:
samples = [sam1, sam2, #sam3, sam4,
            sam5, #sam6,
            sam7, sam8, sam9, sam10,
            sam11, sam12, sam13, sam14, sam15,
            sam16, sam17, sam18, sam19, sam20,
            sam21, sam22, sam23, sam24, sam25,
            sam26, sam27, sam28, sam29, sam30,
            sam31, sam32, sam33, sam34, sam35]

# Quick pre-filtering, these should be low, otherwise it can mess up downstream analysis, but also can get rid of trivial uninteresting things
[sc.pp.filter_cells(sam, min_genes=2000) for sam in samples]
[sc.pp.filter_cells(sam, max_counts=100000) for sam in samples]
[sc.pp.filter_cells(sam, min_counts=5000) for sam in samples]
# Do not filter gene here; concatenate joins on the union, so if a gene fails in a single sample, it will also be deleted from all other samples;

print('Concatenating')
adata = sam1.concatenate(samples[1:])

del samples

adata.X = adata.X.astype('float32')

print(adata)

sc.pl.violin(adata, ['n_genes', 'n_counts'], groupby='replicate', size=0, log=False, cut=0, show=False, save='qc1-pre-norm-replicates.pdf')

# Base filtering for trivial QC failures:
sc.pp.filter_cells(adata, min_genes=3000)
sc.pp.filter_cells(adata, min_counts=8000)
sc.pp.filter_cells(adata, max_counts=100000)
sc.pp.filter_genes(adata, min_cells=50) # Only filter genes here;

print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

#sc.pl.violin(adata, ['n_genes','n_counts'], groupby='stage', size=0, log=False, cut=0, show=False, save='qc1.pdf')
sc.pl.violin(adata, ['n_genes','n_counts'], groupby='replicate', size=0, log=False, cut=0, show=False, save='qc1-replicates.pdf')

p = sb.distplot(adata.obs['n_counts'], kde=False)
p.get_figure().savefig('figures/distplot_ncounts1.pdf')
p = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']<4000], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ncounts2.pdf')
p = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']>10000], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ncounts3.pdf')
#Thresholding decision: genes
p = sb.distplot(adata.obs['n_genes'], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ngenes1.pdf')
p = sb.distplot(adata.obs['n_genes'][adata.obs['n_genes']<2000], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ngenes2.pdf')

print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))

adata.write('./raw_data.h5ad')
