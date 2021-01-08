import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from rpy2.robjects.packages import importr
#from gprofiler import gprofiler
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 1
sc.set_figure_params(dpi=200, dpi_save=300)

sc.settings.figdir = 'markers-leiden0.2'

adata = sc.read('learned.h5ad') #
#sc.pp.log1p(adata)

print(adata.var_names)

oh = open('gene_names.all.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()

marker_genes_dict = {
    'Epiblast': ["Pou5f1"], # Done
    'Primitive Streak': ['Mixl1'], # Done
    'Meso/endoderm': ['Eomes', 'T'], # Done
    'Endoderm': ['Sox17'], # Done
    'Mesoderm': ['Tbx6'], # Done
    'Ectoderm': ['Nr2f1', 'Pax6'],
    'Exe. endoderm': ["Apoa2"], # Done
    'Exe. ectoderm': ["Tfap2c"], # Done
    'Mesenchyme': ['Pmp22'], # Done
    'Blood progenitors': ['Runx1'], # Done
    'Erythroid': ['Gata1'], # Done
    }

sc.pl.stacked_violin(adata, marker_genes_dict, groupby='leiden_r0.2', vmax=3, rotation=90, dendrogram=False, show=False, save='markers.pdf')
sc.pl.dotplot(adata, marker_genes_dict, groupby='leiden_r0.2', dot_max=0.5, dendrogram=False, standard_scale='var', show=False, save='markers.pdf')
sc.pl.heatmap(adata, marker_genes_dict, groupby='leiden_r0.2', vmax=3, show=False, save='markers.pdf')
'''
for k in marker_genes_dict:
    sc.pl.tsne(adata, color=marker_genes_dict[k], size=10, legend_loc='on data', vmax=3, show=False, save='markers-{0}.pdf'.format(k))
    sc.pl.umap(adata, color=marker_genes_dict[k], color_map='plasma', size=10, vmax=3, legend_loc='on data', show=False, save='markers-{0}.pdf'.format(k))

'''
