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

sc.settings.figdir = 'genes'

adata = sc.read('./learned.h5ad')
print(adata)
all_genes = adata.var['n_cells'].index # gene names are stored in the index

TEs = genelist(filename='../../TE_genes_id.mm10.txt', format={'name': 0, 'force_tsv': True})['name']

print(TEs)

for g in all_genes:
    if g not in TEs and '(' not in g:
        print(g)
        sc.pl.umap(adata, color=[g], size=6, legend_loc='on data', color_map='plasma', show=False, save='-{0}.pdf'.format(g), vmin=0, vmax=3)


