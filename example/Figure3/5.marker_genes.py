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

sc.settings.figdir = 'markers'

adata = sc.read('learned.h5ad') # You can skip the script 3 if using te 2b.
#sc.pp.log1p(adata)

print(adata.var_names)

oh = open('gene_names.all.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()

marker_genes_dict = {
    'Epiblast': ["Pou5f1", "Epcam"],
    'Primitive streak': ["Eomes", "Nanog"], #Nanog?!?!
    'Anterior primitive streak': ["Gsc", "Mixl1"],
    'Notochord': ["Noto", "T"],
    'Def. Endoderm': ["Cer1", "Sox7"],
    'Nascent mesoderm': ["Mesp1", "Apela"],
    'Caudal mesoderm': ["Cdx1", "Hes7"],
    'Paraxial mesoderm': ["Tcf15", "Tbx1"],
    'Somitic mesoderm': ["Tbx6", "Dll1"],
    'Pharngyeal mesoderm': ["Tcf21", "Isl1"],
    'Cardiomyocytes': ["Tnnt2", "Myl4"],
    'Allantois': ["Tbx4", "Hoxa11"],
    'Mesenchyme': ["Krt18", "Pmp22"],
    'Hemandothelial prog.': ["Kdr", "Etv2"],
    'Endothelium': ["Pecam1", "Anxa5"],
    'Blood prog.': ["Runx1", "Lmo2"],
    'Erythroid': ["Gata1", "Gypa"],
    'Neuromesoderml prog.': ["Cdx4", "Epha5"],
    'Neurectoderm': ["Six3", "Irx3"],
    'Neural crest': ["Dlx2", "Sox10"],
    'Brain': ["En1", "Pax2"],
    'Spinal cord': ["Sox2", "Pax2"],
    'Surface ectoderm': ["Trp63", "Grhl2"],
    'Visceral endoderm': ["Dkk1", "Amot"],
    'Exe endoderm': ["Ttr", "Apoa2"],
    'Exe ectoderm': ["Tfap2c", "Elf5"],
    'Parietal endoderm': ["Sparc", "Plat"],
    'others': ['Fgf5', 'Lefty2'],
    }

sc.pl.stacked_violin(adata, marker_genes_dict, groupby='leiden_r0.5', rotation=90, dendrogram=True, show=False, save='markers.pdf')
sc.pl.dotplot(adata, marker_genes_dict, groupby='leiden_r0.5', dot_max=0.5, dendrogram=True, standard_scale='var', show=False, save='markers.pdf')
sc.pl.heatmap(adata, marker_genes_dict, groupby='leiden_r0.5', vmax=3, show=False, save='markers.pdf')

for k in marker_genes_dict:
    sc.pl.tsne(adata, color=marker_genes_dict[k], size=10, legend_loc='on data', vmax=3, show=False, save='markers-{0}.pdf'.format(k))
    sc.pl.umap(adata, color=marker_genes_dict[k], color_map='plasma', size=10, vmax=3, legend_loc='on data', show=False, save='markers-{0}.pdf'.format(k))

