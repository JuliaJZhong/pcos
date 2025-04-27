import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import gmean
from pacmap import PaCMAP
from parafac2.normalize import prepare_dataset
from parafac2.parafac2 import parafac2_nd, store_pf2
from sklearn.linear_model import LinearRegression, LogisticRegressionCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

data_path = '/home/jjzhong/projects/pcos/tensor_decomp/data'

def run_pf2(
    X: ad.AnnData,
    rank: int,
    random_state=1,
    doEmbedding: bool = True,
    tolerance=1e-9,
    max_iter: int = 500,
):
    """Run Pf2 model and return modified anndata object"""
    pf_out, _ = parafac2_nd(
        X, rank=rank, random_state=random_state, tol=tolerance, n_iter_max=max_iter
    )

    X = store_pf2(X, pf_out)

    if doEmbedding:
        pcm = PaCMAP(random_state=random_state)
        X.obsm["X_pf2_PaCMAP"] = pcm.fit_transform(X.obsm["projections"])  # type: ignore

    return X

# from https://github.com/meyer-lab/RISE/blob/main/RISE/factorization.py
def correct_conditions(X: ad.AnnData):
    """Correct the conditions factors by overall read depth."""
    sgIndex = X.obs["condition_unique_idxs"]

    counts = np.zeros((np.amax(sgIndex.to_numpy()) + 1, 1))

    cond_mean = gmean(X.uns["Pf2_A"], axis=1)

    x_count = X.X.sum(axis=1)

    for ii in range(counts.size):
        counts[ii] = np.sum(x_count[X.obs["condition_unique_idxs"] == ii])

    lr = LinearRegression()
    lr.fit(counts, cond_mean.reshape(-1, 1))

    counts_correct = lr.predict(counts)

    return X.uns["Pf2_A"] / counts_correct

def logregcv(scoring):
    """Standardizing LogReg for all functions"""
    lrcv = LogisticRegressionCV(
        random_state=0,
        max_iter=10000,
        penalty="l1",   # interesting. also consider elasticnet
        solver="saga",  # faster for larger datasets
        scoring=scoring,
    )
    return lrcv

def run_logregcv(
    adata: ad.AnnData,
    rank: int,
    error_metric: str = "roc_auc",
    max_iter: int = 50,     # for parafac2 decomposition
    test_size: float = 0.2,
    random_state: int = 0
    ):
    
    pf2_output = run_pf2(adata, rank=rank, max_iter=max_iter, doEmbedding=False)
    
    # correct condition factors by read depth
    condition_factors = np.array(correct_conditions(pf2_output))

    # extract disease_status labels
    disease_status = adata.obs.groupby('sample', observed=True)['disease_status'].first()
    labels = disease_status.values

    # convert from string to binary values
    label_mapping = {"HC": 0, "PCOS": 1}
    labels = np.array([label_mapping[label] for label in labels])

    X_train, X_test, y_train, y_test = train_test_split(
        condition_factors, 
        labels, 
        test_size=test_size, 
        random_state=random_state, 
        stratify=labels)
    
    print('y_train:', y_train)
    print('y_test:', y_test)

    clf = logregcv(error_metric).fit(X_train, y_train)

    if error_metric == "roc_auc":
        confidence_scores = clf.decision_function(X_test)
        print('confidence scores:', confidence_scores)
        score = roc_auc_score(y_test, confidence_scores)

    # elif error_metric == "accuracy": TODO
        score = clf.score(X_test, y_test)
        print('accuracy:', score)

    return score

adata = sc.read_h5ad(data_path + '/integrated_2025-04-26.h5ad')
adata = prepare_dataset(adata, condition_name='sample', geneThreshold=0.01) 
adata.obs["condition_unique_idxs"] = adata.obs['condition_unique_idxs'].astype('category')

ranks_to_test = [10, 30] # [10, 30, 50, 80, 100]  
score_results = np.zeros(len(ranks_to_test))

for idx, rank in enumerate(ranks_to_test):
    print(f'\n\n # components: {rank}')
    score_results[idx] = run_logregcv(adata, rank=rank, error_metric='roc_auc')

print('ranks tested:', ranks_to_test)
print('roc auc scores:', score_results)
