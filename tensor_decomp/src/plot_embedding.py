import anndata as ad
import pandas as pd
import seaborn as sns
import scanpy as sc

data_path = '/home/jjzhong/projects/pcos/tensor_decomp/data'
fig_path = '/home/jjzhong/projects/pcos/tensor_decomp/figures'

cmap = sns.diverging_palette(240, 10, center="dark", as_cmap=True)

pf2 = sc.read_h5ad(data_path + '/pf2_stacas_umap_2025-04-26.h5ad')

# add cell annotations to pf2 anndata object
cell_annot = pd.read_csv(data_path + '/cell_annot.csv')
cell_annot = cell_annot.set_index("barcode")
shared_cells = cell_annot.loc[cell_annot.index.intersection(pf2.obs_names)]
pf2.obs['cell_annot'] = pf2.obs_names.map(shared_cells['cell_type'])

sc.pp.neighbors(pf2, use_rep='X_pf2_UMAP')

# initial UMAP with color -> cell type
sc.pl.umap(
    pf2,
    color='cell_annot',
    # setting a smaller point size to get prevent overlap
    size=2,
    # legend_loc="on data"
    save='_pf2_cell_annot.png'
)

# UMAPs colored by each component's weighted projections
wp = pf2.obsm['weighted_projections']

for comp in range(wp.shape[1]):
    pf2.obs['WP_' + str(comp)] = wp[:, comp]

    sc.pl.umap(
    pf2,
    color='WP_' + str(comp),
    color_map=cmap,
    vmin=-0.2,
    vmax=0.2,
    vcenter=0,
    # setting a smaller point size to get prevent overlap
    size=2,
    # legend_loc="on data"
    save='_pf2_WP_' + str(comp) + '.png'
    )