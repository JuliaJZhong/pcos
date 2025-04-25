import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scanpy as sc

data_path = '/home/jjzhong/projects/pcos/tensor_decomp/data'
fig_path = '/home/jjzhong/projects/pcos/tensor_decomp/figures'

cmap = sns.diverging_palette(240, 10, as_cmap=True)

adata = sc.read_h5ad(data_path + '/pf2_2025-04-22.h5ad')

# note: weighted_projections = projections x eigen-state factors
component_weights = adata.uns['Pf2_weights']
projections = adata.obsm['projections']
weighted_projections = adata.obsm['weighted_projections']

condition_factors = adata.uns['Pf2_A']
eigen_factors = adata.uns['Pf2_B']
gene_factors = adata.varm['Pf2_C']

# plotting histograms of distributions
hists = plt.figure(figsize=(25, 8))
ax1 = hists.add_subplot(1, 3, 1)
conditions, rank = condition_factors.shape
sns.histplot(
condition_factors, 
    ax=ax1, 
)
ax1.set_title('Condition factor distribution')

# eigen-state factors
ax2 = hists.add_subplot(1, 3, 2)
eigen, rank = eigen_factors.shape
sns.histplot(
    eigen_factors, 
    ax=ax2, 
    )
ax2.set_title('Eigen-state factor distribution')

# gene factors
ax3 = hists.add_subplot(1, 3, 3)
gene, rank = gene_factors.shape
sns.histplot(
    gene_factors, 
    ax=ax3, 
    )
ax3.set_title('Gene factor distribution')

hists.savefig(fig_path + "/histograms.png")

# --------
# plotting heatmaps
# conditions factors
# TODO: add sample labels. these may not be correct. also add that colorbar thing. also hierarchical
heatmaps = plt.figure(figsize=(25, 8))
ax1 = heatmaps.add_subplot(1, 3, 1)
conditions, rank = condition_factors.shape
sns.heatmap(
    condition_factors, 
    ax=ax1, 
    xticklabels=range(1, rank + 1), 
    yticklabels=list(adata.obs['sample'].cat.categories),
    cmap=cmap
)
ax1.set_xlabel('Components')
ax1.set_ylabel('Samples')
ax1.set_title('Condition factors')

# eigen-state factors
ax2 = heatmaps.add_subplot(1, 3, 2)
eigen, rank = eigen_factors.shape
sns.heatmap(
    eigen_factors, 
    ax=ax2, 
    xticklabels=range(1, rank + 1), 
    yticklabels=range(1, eigen + 1),
    cmap=cmap
    )
ax2.set_xlabel('Components')
ax2.set_ylabel('Eigen-state')
ax2.set_title('Eigen-state factors')

# gene factors
# TODO: need to filter these
# TODO: add gene labels
# TODO: the loadings are so small. check distributions of these? ok i think we fine actually
# TODO: hierarchical
ax3 = heatmaps.add_subplot(1, 3, 3)
gene, rank = gene_factors.shape
sns.heatmap(
    gene_factors, 
    ax=ax3, 
    xticklabels=range(1, rank + 1), 
    yticklabels=range(1, gene + 1), 
    cmap=cmap
    )
ax3.set_xlabel('Components')
ax3.set_ylabel('Gene')
ax3.set_title('Gene factors')

heatmaps.savefig(fig_path + "/heatmaps.png")