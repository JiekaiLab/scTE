import logging, matplotlib, os, sys
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors

from glbase3 import *

plt.rcParams['figure.figsize']=(8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=10

sc.settings.figdir = 'tes'

adata = sc.read('./learned.h5ad')
print(adata)
all_genes = adata.var['n_cells'].index # gene names are stored in the index

TEs = genelist(filename='TE_genes_id.mm10.txt.gz', format={'name': 0, 'force_tsv': True}, gzip=True)

#merker_tes = ['ID2', 'MER5C1', 'MER34B-int', 'MER63D', 'MT2A']
#sc.pl.stacked_violin(adata, var_names=merker_tes, groupby='leiden_r0.2', rotation=90, show=False, save='tes.pdf')

for te in TEs:
    print(te['name'])
    if te['name'] in all_genes:
        sc.pl.umap(adata, color=[te['name'], te['name']], size=10, legend_loc='on data', show=False, save='TE-{0}.pdf'.format(te['name']), vmin=0, vmax=3)


