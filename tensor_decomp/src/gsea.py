# base code by favour
import anndata as ad
import pandas as pd
import scanpy as sc
import gseapy as gp
import numpy as np
from gseapy.plot import gseaplot

def parse_gmt(file_path):
    gene_sets = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            gene_set_name = parts[0]
            description = parts[1]
            genes = parts[2:]
            gene_sets[gene_set_name] = genes
    return gene_sets

"""DGE Analysis b/w conditions"""

data_path = '/home/jjzhong/projects/pcos/tensor_decomp/data'
fig_path = '/home/jjzhong/projects/pcos/tensor_decomp/figures'

adata = sc.read_h5ad(data_path + '/pf2_stacas_2025-04-26.h5ad')
gene_factors = np.array(adata.varm['Pf2_C'])

# look specifically at component with highest weight first
component = np.argmax(adata.uns['Pf2_weights'])

# list of genes and their corresponding factor values for above selected component
ranked_genes = pd.DataFrame()
ranked_genes['gene_symbols'] = adata.var_names.values
ranked_genes['latent_factors'] = gene_factors[:, component]

# path to the full C2 GMT file you downloaded
c2_gmt_file = data_path + '/c2.all.v2024.1.Hs.symbols.gmt'

# load the gene sets
all_gene_sets = parse_gmt(c2_gmt_file)

# filter gene sets with keywords like 'steroid', 'cholesterol'
keywords = ['steroid', 'cholesterol', 'lipid', 'inflamm', 'insulin', 'ovar', 'androgen']
filtered_gene_sets = {name: genes for name, genes in all_gene_sets.items()
                      if any(k in name.lower() for k in keywords)}

print(f"Found {len(filtered_gene_sets)} matching gene sets:")
for name in filtered_gene_sets:
    print(name)

# run preranked GSEA
pre_res = gp.prerank(
    rnk=ranked_genes,                                   # dataframe of genes and scores    
    gene_sets=filtered_gene_sets,                       # your custom gene set file
    outdir= fig_path + '/GSEA',                         # output folder
    min_size=5,
    max_size=500,
    permutation_num=1000,                              
    seed=6,
    threads=4,
    verbose=True
)