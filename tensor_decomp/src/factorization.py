import anndata
import cupy as cp
from datetime import date
import numpy as np
import pandas as pd
from pacmap import PaCMAP   # TODO: remove
import parafac2 as pf2
from parafac2.normalize import prepare_dataset
from parafac2.parafac2 import parafac2_nd, store_pf2
import scanpy as sc
from sklearn.preprocessing import StandardScaler
import umap

data_path = '/home/jjzhong/projects/pcos/tensor_decomp/data'

print('number of CUDA devices:', cp.cuda.runtime.getDeviceCount())
print('-----')

# tweak from https://github.com/meyer-lab/RISE/blob/main/RISE/factorization.py
def run_pf2(
    X: anndata.AnnData,
    rank: int,
    random_state=1,
    doEmbedding: bool = True,
    tolerance=1e-9,
    max_iter: int = 500,
):
    """Run Pf2 model and store results in anndata file"""
    pf_out, _ = parafac2_nd(
        X, rank=rank, random_state=random_state, tol=tolerance, n_iter_max=max_iter
    )

    X = store_pf2(X, pf_out)

    if doEmbedding:
        scaled = StandardScaler().fit_transform(X.obsm["projections"])
        reducer = umap.UMAP()
        X.obsm["X_pf2_UMAP"] = reducer.fit_transform(scaled)

    return X

# load the filtered and normalized output of preprocessing.py
adata = sc.read_h5ad(data_path + '/integrated_2025-04-26.h5ad')

adata = prepare_dataset(adata, condition_name='sample', geneThreshold=0) # TODO: 0.01
print('adata object:', adata)
print('-----')
adata.obs["condition_unique_idxs"] = adata.obs['condition_unique_idxs'].astype('category')
print('adata.obs["condition_unique_idxs"]:', adata.obs['condition_unique_idxs'])
print('-----')
print('adata.obs_vector("sample"):', adata.obs_vector('sample'))
print('-----')

adata = run_pf2(adata, rank=30, max_iter=50)

print('integrated adata object post-pf2:', adata)

today = date.today()
adata.write_h5ad(filename=data_path + '/pf2_stacas_' + str(today) + '.h5ad' )