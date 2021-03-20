import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
#from rpy2.robjects.packages import importr
#from gprofiler import gprofiler
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 1
sc.set_figure_params(dpi=200, dpi_save=300)

#matplotlib.rcParams['pdf.fonttype']=42
#matplotlib.rcParams['font.size']=6

todo = 'leiden_r0.3'

sc.settings.figdir = 'markers-{0}'.format(todo)

adata = sc.read('learned.h5ad')

marker_genes_dict = {
    'Epiblast': ["Pou5f1"],
    'Primitive streak': ["Mixl1"], #Nanong?!?!
    'Endoderms': ["Cer1", "Sox7"],
    'Mesoderms': ["T", 'Cdx1'],
    'Ectoderms': ['Six3'], # And Grhl2

    'Exe endoderm': ["Apoa2"],
    'Exe ectoderm': ["Tfap2c"],

    'Cardiomyocytes': ["Tnnt2"],
    'Blood prog.': ["Lmo2", ],
    'Erythroid': ["Gypa"],
    }

sc.pl.stacked_violin(adata, marker_genes_dict, groupby=todo, rotation=90, dendrogram=True, show=False, save='markers.pdf')
sc.pl.dotplot(adata, marker_genes_dict, groupby=todo, color_map='Greens', dot_max=0.7, dendrogram=True, standard_scale='var', show=False, save='markers.pdf')
sc.pl.heatmap(adata, marker_genes_dict, groupby=todo, vmax=3, show=False, save='markers.pdf')

# high, few: Expressed rarely, but very high in the cells that they are expressed in
marker_genes_dictB = {
    #'Epiblast': ['MTEb-int',],
    #'Primitive streak': ['RLTR1D2_MM', ],
    #'Endothelium': ['ERVB7_2B-LTR_MM',],

    #'Ectoderms': ['MamRep137'],
    #'Endoderms': ['MLT1I'],
    'Mesoendoderm': ['RLTR48A', 'IAPEY4_LTR', 'ORR1F-int'],
    'Extraembryonic': ['LTR16A', ],
    'Exe. endoderm': ['MER5C', 'RLTR6B_Mm',],
    #'Exe. ectoderm': ['ERVB4_2-LTR_MM', ],
    'Cardiomyocyte': ['L1ME3D', 'RLTR13A2', 'ERVB2_1A-I_MM-int', 'RLTR16'],
    }
sc.pl.dotplot(adata, marker_genes_dictB, groupby=todo, dot_max=0.3, dendrogram=True, standard_scale='var', vmax=1, show=False, save='markersB.pdf')

# Super-specific
marker_genes_dictC = {
    #'Primitive streak': [ ],
    'Mesoendoderm': ['ERVB4_1C-LTR_Mm', 'ETnERV3-int',],
    #'others':['MuRRS4-int'],
    'Exe. endoderm': ['MER46C', 'MuRRS4-int',  'RLTR20B3',  'RLTR1B-int', 'LTRIS2',],
    'Exe. ectoderm': ['RLTR45', 'RLTR45-int', 'IAPLTR1_Mm'],
    #'Cardiomyocyte': ['ETnERV3-int', 'L1ME3D', 'RLTR13A2', 'ERVB2_1A-I_MM-int'],
    'Erythroid': ['RLTR10F', 'L1_Mur1',],
    }
sc.pl.dotplot(adata, marker_genes_dictC, groupby=todo, dot_max=0.7, dendrogram=True, standard_scale='var', vmax=1, show=False, save='markersC.pdf')
