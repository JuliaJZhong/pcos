import anndata
import cupy as cp
from datetime import date
import numpy as np
import pandas as pd
from pacmap import PaCMAP
import parafac2 as pf2
from parafac2.normalize import prepare_dataset
from parafac2.parafac2 import parafac2_nd, store_pf2
import scanpy as sc

data_path = '/home/jjzhong/tensor_decomp/data'

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
        pcm = PaCMAP(random_state=random_state)
        X.obsm["X_pf2_PaCMAP"] = pcm.fit_transform(X.obsm["projections"])  # type: ignore

    return X

# load the filtered and normalized output of preprocessing.py
adata = sc.read_h5ad(data_path + '/adata_2025-04-20.h5ad')

# TODO: deeper dive into geneThreshold
adata = prepare_dataset(adata, condition_name='sample', geneThreshold=0.01)  
print('adata object:', adata)
print('-----')
adata.obs["condition_unique_idxs"] = adata.obs['condition_unique_idxs'].astype('category')
print('adata.obs["condition_unique_idxs"]:', adata.obs['condition_unique_idxs'])
print('-----')
print('adata.obs_vector("sample"):', adata.obs_vector('sample'))
print('-----')

# TODO run pf2 function. pick rank for it. and plot pacmap. also pull and push from github. might need to make new keys
adata = run_pf2(adata, rank=30, max_iter=100)

print('adata object post-pf2:', adata)

today = date.today()
adata.write_h5ad(filename=data_path + '/pf2_' + str(today) + '.h5ad' )