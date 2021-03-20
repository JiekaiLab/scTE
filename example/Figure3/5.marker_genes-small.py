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

sc.settings.figdir = 'markers-small'

adata = sc.read('learned.h5ad')

marker_genes_dict = {
    'Epiblast': ["Pou5f1"],
    'Primitive streak': ["Eomes", "Mixl1"], #Nanong?!?!
    'Endoderms': ["Cer1", "Sox7"],
    'Mesoderms': ["T", 'Cdx1'],
    'Ectoderms': ['Grhl2', 'Six3'],

    'Exe endoderm': ["Apoa2"],
    'Exe ectoderm': ["Tfap2c"],

    'Cardiomyocytes': ["Tnnt2"],
    'Blood prog.': ["Lmo2", ],
    'Erythroid': ["Gypa"],
    }

sc.pl.stacked_violin(adata, marker_genes_dict, groupby='leiden_r0.5', rotation=90, dendrogram=True, show=False, save='markers.pdf')
sc.pl.dotplot(adata, marker_genes_dict, groupby='leiden_r0.5', color_map='Greens', dot_max=0.5, dendrogram=True, standard_scale='var', show=False, save='markers.pdf')
sc.pl.heatmap(adata, marker_genes_dict, groupby='leiden_r0.5', vmax=3, show=False, save='markers.pdf')

for k in marker_genes_dict:
    sc.pl.tsne(adata, color=marker_genes_dict[k], size=10, legend_loc='on data', vmax=3, show=False, save='markers-{0}.pdf'.format(k))
    sc.pl.umap(adata, color=marker_genes_dict[k], color_map='plasma', size=10, vmax=3, legend_loc='on data', show=False, save='markers-{0}.pdf'.format(k))
