import logging, matplotlib, os, sys
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import pandas as pd
from glbase3 import genelist
plt.rcParams['figure.figsize']=(8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=10

sc.settings.figdir = 'diffexp'

adata = sc.read('./learned.h5ad')

sc.tl.rank_genes_groups(adata, 'leiden_r0.5', method='wilcoxon', n_genes=3000)
adata.write('./de.h5ad')

adata = sc.read('./de.h5ad')

sc.pl.rank_genes_groups(adata, n_genes=25, sharey=True, show=False, save='genes-top25.pdf')
sc.pl.rank_genes_groups(adata, key='rank_genes_groups', show=False, save='genes.pdf')
sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups', show=False, save='genes-top25.pdf')

#print(pd.DataFrame(adata.uns['rank_genes_groups']))

print(pd.DataFrame(adata.uns['rank_genes_groups']['names']))

print()
topall = pd.DataFrame(adata.uns['rank_genes_groups']['names']) # get all;
fcs = pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges'])
padj = pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj'])

topall.to_csv('top100.csv')

# Go through and trim the TEs:

TEs = set(genelist(filename='../../TE_genes_id.mm10.txt', format={'name': 0, 'force_tsv': True})['name'])

newcols = {}

groups = list(topall.columns.values)

for group in groups:
    newcols[group] = []

    t = zip([i[group] for i in adata.uns['rank_genes_groups']['names']], [i[group] for i in adata.uns['rank_genes_groups']['logfoldchanges']], [i[group] for i in adata.uns['rank_genes_groups']['pvals_adj']])

    print('Group: {0}'.format(group))
    print(t)

    for item in t:
        print(item)
        if abs(item[1]) < 1: # fold change
            continue
        if item[2] > 0.01: # just in case
            continue

        if item[0] in TEs:
            newcols[group].append(item[0])


# join all and draw a dotplot:
joined = []
for group in newcols:
        joined += newcols[group]

# Need to remove duplicates, but preserver order:
newl = []
for i in joined:
    if i not in newl:
        newl.append(i)
joined = newl

print(joined)
sc.pl.dotplot(adata, joined, groupby='leiden_r0.5', dot_max=0.7, dendrogram=True, standard_scale='var', show=False, save='de-tes.pdf')
sc.pl.matrixplot(adata, joined, groupby='leiden_r0.5', dendrogram=True, standard_scale='var', show=False, save='de-tes.pdf')

for k in joined:
    sc.pl.tsne(adata, color=[k,k], size=15, legend_loc='on data', vmax=2, show=False, save='markers-{0}.pdf'.format(k))
    sc.pl.umap(adata, color=[k,k], size=15, legend_loc='on data', vmax=2, show=False, save='markers-{0}.pdf'.format(k))
