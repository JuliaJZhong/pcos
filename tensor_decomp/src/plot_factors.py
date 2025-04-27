import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import seaborn as sns
import scanpy as sc
from matplotlib.axes import Axes
from matplotlib.patches import Patch

data_path = '/home/jjzhong/projects/pcos/tensor_decomp/data'
fig_path = '/home/jjzhong/projects/pcos/tensor_decomp/figures'

cmap = sns.diverging_palette(240, 10, as_cmap=True)

# from https://github.com/meyer-lab/RISE/blob/main/RISE/figures/commonFuncs/plotFactors.py
def reorder_table(projs: np.ndarray):
    '''
    Hierarchical clustering for reordering tables. Used for condition factors and gene factors
    '''
    assert projs.ndim == 2
    linkage_mtx = sch.linkage(projs, method="complete", metric="cosine", optimal_ordering=True)
    return sch.leaves_list(linkage_mtx)

def plot_histograms(adata: ad.AnnData):
    '''
    Plot distributions of:
    - log_10 transformed conditions factors
    - eigen-state factors
    - log_10 transformed gene factors
    '''
    condition_factors = np.array(adata.uns['Pf2_A'])
    eigen_factors = np.array(adata.uns['Pf2_B'])
    gene_factors = np.array(adata.varm['Pf2_C'])

    condition_factors_log10 = np.log10(condition_factors)
    gene_factors_log10 = np.log10(gene_factors)

    # condition factors
    hists = plt.figure(figsize=(25, 8))
    ax1 = hists.add_subplot(1, 3, 1)
    sns.histplot(
        condition_factors_log10.flatten(), 
        ax=ax1, 
    )
    ax1.set_title('Condition factor distribution')
    ax1.set_xlabel(r'$log_{10}$(loading)')

    # eigen-state factors
    ax2 = hists.add_subplot(1, 3, 2)
    sns.histplot(
        eigen_factors.flatten(), 
        ax=ax2, 
    )
    ax2.set_title('Eigen-state factor distribution')
    ax2.set_xlabel('loading')

    # gene factors
    ax3 = hists.add_subplot(1, 3, 3)
    sns.histplot(
        gene_factors_log10.flatten(), 
        ax=ax3, 
    )
    ax3.set_title('Gene factor distribution')
    ax3.set_xlabel(r'$log_{10}$(loading)')

    hists.savefig(fig_path + "/histograms_stacas.png")

def plot_conditions(adata: ad.AnnData, ax: Axes):
    '''
    Standardize, hierarchically cluster, and plot heatmap of condition factor matrix
    '''
    condition_factors = np.array(adata.uns['Pf2_A'])
    conditions, rank = condition_factors.shape
    sample_ids = adata.obs['sample'].cat.categories
    
    # log-transform condition factors
    X = np.log(condition_factors)

    # Z-score down components/columns
    X -= np.mean(X, axis=0)
    X /= np.std(X, axis=0)

    # hierarchical clustering
    idxs = reorder_table(X)
    X = X[idxs]
    sample_ids = sample_ids[idxs]

    # disease status labels
    disease_status = adata.obs.groupby('sample', observed=True)['disease_status'].first()
    labels = disease_status.iloc[idxs].values

    # disease status colors
    colors = sns.color_palette(n_colors=pd.Series(labels).nunique()).as_hex()

    # legend
    group_2_color = {}
    legend_elements = []
    for index, group in enumerate(pd.Series(labels).unique()):
        group_2_color[group] = colors[index]
        legend_elements.append(Patch(color=colors[index], label=group))
    ax.legend(handles=legend_elements, bbox_to_anchor=(0.18, 1.07))

    # drawing disease status bar on y-axis
    row_colors = pd.Series(labels).map(group_2_color)
    for rows, color in enumerate(row_colors):
        ax.add_patch(
            plt.Rectangle(
                xy=(-0.05, rows),
                width=0.05,
                height=1,
                color=color,
                lw=0,
                transform=ax.get_yaxis_transform(),
                clip_on=False,
            )
        )

    sns.heatmap(
        X, 
        ax=ax, 
        xticklabels=range(1, rank + 1), 
        yticklabels=list(sample_ids),
        cmap=cmap,
        cbar_kws={'label': r'Z-score of Log$_{10}$(loading)'}
    )
    ax.yaxis.set_tick_params(pad=25, left=False)
    ax.set_xlabel('Components')
    ax.set_ylabel('Samples')
    ax.set_title('Condition factors')

def plot_eigens(adata: ad.AnnData, ax: Axes):
    '''
    Scale and plot heatmap of eigen-state factor matrix
    '''
    eigen_factors = np.array(adata.uns['Pf2_B'])
    eigen, rank = eigen_factors.shape

    # scale factors to [-1, 1]
    X = eigen_factors / np.max(np.abs(eigen_factors))

    sns.heatmap(
        X, 
        ax=ax, 
        xticklabels=range(1, rank + 1), 
        yticklabels=range(1, eigen + 1),
        cmap=cmap,
        center=0, 
        vmin=-1,
        vmax=1,
        cbar_kws={'label': r'loading/loading$_{max}$'}
    )
    ax.set_xlabel('Components')
    ax.set_ylabel('Eigen-state')
    ax.set_title('Eigen-state factors')

def plot_genes(adata: ad.AnnData, ax: Axes, thresh=0):
    '''
    Filter to top contributing genes, hierarchically cluster, scale, and plot heatmap of gene factor matrix
    '''
    gene_factors = np.array(adata.varm['Pf2_C'])
    gene, rank = gene_factors.shape
    
    X = gene_factors
    gene_symbols = adata.var.index.values

    # filter out to genes who have a contribution > thresh to at least 1 component
    max_weights = np.max(np.abs(X), axis=1)
    keep = max_weights > thresh
    X = X[keep]
    gene_symbols = gene_symbols[keep]

    # hierarchical clustering
    idxs = reorder_table(X)
    X = X[idxs]
    gene_symbols = gene_symbols[idxs]

    # scale factors to [-1, 1]
    X = X / np.max(np.abs(X))
    
    sns.heatmap(
        X, 
        ax=ax, 
        xticklabels=range(1, rank + 1), 
        yticklabels=gene_symbols, 
        cmap=cmap,
        center=0, 
        vmin=-1,
        vmax=1,
        cbar_kws={'label': r'loading/loading$_{max}$'}
    )
    ax.set_xlabel('Components')
    ax.set_ylabel('Gene')
    ax.set_title('Gene factors')

def plot_heatmaps(adata: ad.AnnData, gene_thresh=0):
    ''' 
    Plot heatmaps of:
    - hierarchically-clustered log_10 transformed conditions factors
    - eigen-state factors
    - hierarchically-clustered gene factors (filtered to genes that have factor/loading > gene_thresh for at least 1 component)
    '''
    heatmaps = plt.figure(figsize=(30, 8))

    ax1 = heatmaps.add_subplot(1, 3, 1)
    plot_conditions(adata, ax1)

    ax2 = heatmaps.add_subplot(1, 3, 2)
    plot_eigens(adata, ax2)

    ax3 = heatmaps.add_subplot(1, 3, 3)
    plot_genes(adata, ax3, thresh=gene_thresh)

    heatmaps.savefig(fig_path + "/heatmaps_stacas.png")

adata = sc.read_h5ad(data_path + '/pf2_stacas_2025-04-26.h5ad')
plot_histograms(adata)
plot_heatmaps(adata, gene_thresh=0.13)  # gene_thresh only for visualization purposes