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
matplotlib.rcParams['font.size']=6

sc.settings.figdir = 'specific-tes'

adata = sc.read('./learned.h5ad')

# high, few: Expressed rarely, but very high in the cells that they are expressed in
marker_genes_dictB = {
    #'Epiblast': ['MTEb-int',],
    'Primitive streak': ['RLTR1D2_MM', ],
    #'Endothelium': ['ERVB7_2B-LTR_MM',],

    #'Ectoderms': ['MamRep137'],
    #'Endoderms': ['MLT1I'],
    'Mesoendoderm': ['RLTR48A', 'IAPEY4_LTR', 'ORR1F-int'],
    'Extraembryonic': ['LTR16A', ],
    'Exe. endoderm': ['MER5C', 'RLTR6B_Mm',],
    #'Exe. ectoderm': ['ERVB4_2-LTR_MM', ],
    'Cardiomyocyte': ['L1ME3D', 'RLTR13A2', 'ERVB2_1A-I_MM-int', 'RLTR16'],
    }
sc.pl.dotplot(adata, marker_genes_dictB, groupby='leiden_r0.5', dot_max=0.3, dendrogram=True, standard_scale='var', vmax=1, show=False, save='markersB.pdf')

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
sc.pl.dotplot(adata, marker_genes_dictC, groupby='leiden_r0.5', dot_max=0.7, dendrogram=True, standard_scale='var', vmax=1, show=False, save='markersC.pdf')
