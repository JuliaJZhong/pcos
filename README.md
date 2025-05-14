# Single-cell transcriptomics analysis for elucidation of molecular signatures of polycystic ovary syndrome ⚓️ 👩‍👩‍👧‍👦 🤺 🍄 🧬
Polycystic ovary syndrome (PCOS) is an under-researched hormonal disease that affects an estimated 6–13% of women of reproductive age worldwide [1]. Understanding the molecular mechanisms underlying the pathophysiology of the disease is crucial for developing novel diagnostics and targeted treatments. In this work, we analyzed a single-cell RNA sequencing dataset of ovarian theca cells obtained from healthy controls and PCOS patients via clustering and cell-type annotation, tensor decomposition via Parallel Factor Analysis 2 (PARAFAC2), and gene set enrichment analysis (GSEA) to elucidate upregulation and downregulation of enriched pathways. Building upon the original authors’ analysis with computational methods that harness cell-type specificity, we pinpointed potential areas of dysregulation in PCOS surrounding ovarian insufficiency, inflammatory and oxidative stress response, lipid metabolism, and insulin growth factor signaling. 

# Overview 🌈
This repository contains code for analyzing scRNA-seq data of cells extracted from the ovarian theca interna of PCOS patients and healthy controls. Analysis includes preprocessing (including Sub-Type Anchor Correction for Alignment in Seurat (STACAS)), clustering and cell-type annotation, PARAFAC2 tensor decomposition, and gene set enrichment analysis (GSEA).

Python and R ibraries used for analysis are cited at the bottom of this page.

# Data 📊
The dataset we analyzed consists of single-cell RNA sequencing data collected from in vitro culture of ovarian theca cells from 10 patients (5 PCOS, 5 healthy controls). "Human theca interna tissue was obtained from follicles of women undergoing hysterectomy". After treatment, processing, and freezing in liquid N2, all samples were sent to "Active Motif (Carlsbad, CA, USA) for... single-cell library preparation, using Active Motifs proprietary conditions. Following single-cell library preparation, 10× single-cell RNA sequencing (scRNA-seq) was performed using an Illumina NextSeq 500 (San Diego, CA, USA) sequencing apparatus to generate 91 bp sequencing reads"[2].

You can access the .zip file of raw data online [here](https://zenodo.org/records/7942968). You may contact us for direct access to the .h5ad file that contains the integrated, preprocessed data stored as an `anndata` object (you may also run the code in the preprocessing scripts). 

# Structure 🌲
```
annotation/
└── clustering_and_annotation.py       # script for Leiden clustering and cell-type annotation

gsea/
├── gsea.py                            # script for GSEA between PCOS and healthy control groups   
└── pf2_gsea.py                        # script for GSEA on PARAFAC2 gene factors

preprocessing/
├── 10x_to_anndata.py                  # converting feature-barcode matrices into an anndata object
├── data_integration/
│   ├── pcos.Rproj                     # R project configuration file
│   ├── renv                           
│   │   ├── activate.R                 # script to activate R environment
│   │   └── settings.json              # stores Quarto project metadata
│   ├── renv.lock                      # R project package information
│   └── stacas.qmd                     # Quarto markdown notebook for STACAS data integration
├── preprocessing.py                   # script for quality control, normalization, feature selection, dimensionality reduction
└── unzip.py                           # script for unzipping original .zip file of feature-barcode matrices

tensor_decomp/
├── info.txt                           # miscellaneous notes regarding data files
└── src
    ├── factorization.py               # script for performing PARAFAC2 tensor decomposition
    ├── log_reg.py                     # script for performing logistic regression CV with LASSO regularization on a range of PARAFAC2 decompositions
    ├─- plot_embedding.py              # script for plotting UMAP embeddings overlayed with PARAFAC2 weighted projections
    ├── plot_factors.py                # script for plotting heatmaps of sample factors, eigen-state factors, gene factors
    └── plot_triangle.py               # script for performing un-penalized logistic regression on pairs of components and plotting heatmap of prediction accuracy
```

# Installation 🔧
To run the code, first, download the script(s) of your choice. Make sure the correct versions of each package are installed in your environment. Be sure to change `data_path`and `fig_path` variables to your data and desired figure output paths respectively. Now, you may run the script. 

Note that to install the Python package [parafac2](https://github.com/meyer-lab/parafac2/tree/main#) from GitHub, run the following in your virtual environment:
```
pip install git+https://github.com/meyer-lab/parafac2.git@main
```
To install the R packages [schard](https://github.com/cellgeni/schard/tree/main) and [STACAS](https://github.com/carmonalab/STACAS) from GitHub, you can use [remotes](https://github.com/r-lib/remotes):
```
install.packages("remotes")
remotes::install_github("cellgeni/schard")
remotes::install_github("carmonalab/STACAS")
```

Software & Packages
```
anndata==0.11.4
cupy==13.4.1
gseapy==1.1.8
h5py==3.7
pacmap==0.8.0
parafac2==1.0.0
scanpy==1.11.1
scikit_learn==1.6.1
scipy==1.15.3
seaborn==0.13.2
umap_learn==0.5.7
```

For gsea.py and clustering_and_annotation.py:
```
python==3.9
matplotlib==3.10.0
numpy==2.0.2
pandas==2.2.2
```
For preprocessing/ Python scripts, tensor_decomp/, and pf2_gsea.py:
```
python==3.13.3
matplotlib==3.10.3
numpy==2.2.5
pandas==2.2.3
```
For preprocessing/data_integration/:
```
anndata==0.7.5.6
dplyr==1.1.4
ggplot2==3.5.2
schard==0.0.1
Seurat==5.3.0
STACAS==2.2.2                                       
```

# References 📋
[1] World Health Organization. Polycystic Ovary Syndrome. WHO, 2025.

[2] Harris RA, McAllister JM, Strauss JF 3rd. Single-Cell RNA-Seq Identifies Pathways and Genes Contributing to the Hyperandrogenemia Associated with Polycystic Ovary Syndrome. Int J Mol Sci. 2023 Jun 25;24(13):10611. doi: 10.3390/ijms241310611. PMID: 37445796; PMCID: PMC10341507.
