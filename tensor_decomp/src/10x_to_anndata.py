import anndata as ad
import os
import scanpy as sc

data_path = '/home/jjzhong/projects/pcos/tensor_decomp/data'
raw_data_path = os.path.join(data_path, 'raw_data')

samples = {} # dictionary with key = sample IDs, value = data path

for sample_id in os.listdir(raw_data_path):
  sample_path = os.path.join(raw_data_path, sample_id)
  subdir = os.listdir(sample_path)                              # i know there is only 1 subdir
  subsubdir = os.listdir(os.path.join(sample_path, subdir[0]))  # i know there is only 1 subsubdir

  path = os.path.join(sample_path, subdir[0], subsubdir[0])
  samples[sample_id] = path

adatas = {} # dictionary with key = sample IDs, value = AnnData object of scRNAseq data

for sample_id, path in samples.items():
  adata = sc.read_10x_mtx(path, var_names='gene_symbols', make_unique=True, cache=True)
  adatas[sample_id] = adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
print('sample value counts:', adata.obs["sample"].value_counts())
print('-----')
print('initial anndata object:', adata)

adata.write_h5ad(filename=data_path + '/raw_adata.h5ad' )