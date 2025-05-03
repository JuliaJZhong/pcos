import anndata as ad
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

data_path = '/home/jjzhong/projects/pcos/tensor_decomp/data'
fig_path = '/home/jjzhong/projects/pcos/tensor_decomp/figures'

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

pf2_output = sc.read_h5ad(data_path + '/pf2_stacas_umap_2025-04-26.h5ad')

# correct condition factors by read depth
condition_factors = np.array(correct_conditions(pf2_output))

# extract disease_status labels
disease_status = pf2_output.obs.groupby('sample', observed=True)['disease_status'].first()
labels = disease_status.values

# convert from string to binary values
label_mapping = {"HC": 0, "PCOS": 1}
labels = np.array([label_mapping[label] for label in labels])

num_components = condition_factors.shape[1]

# logistic regression scores matrices
roc_aucs = np.zeros((num_components, num_components))
accuracies = np.zeros((num_components, num_components))

# evaluate performance of pairs of components in predicting disease status
for c in range(num_components):
    for r in range(c, num_components):

        # concatenate the two components of interest
        data = condition_factors[:, [r, c]]

        X_train, X_test, y_train, y_test = train_test_split(
        data, 
        labels, 
        test_size=0.4, 
        random_state=0, 
        stratify=labels)

        print('y_train:', y_train)
        print('y_test:', y_test)

        clf = LogisticRegression(penalty=None, random_state=0).fit(X_train, y_train)

        confidence_scores = clf.decision_function(X_test)
        roc_aucs[r, c] = roc_auc_score(y_test, confidence_scores)

        accuracies[r, c] = clf.score(X_test, y_test)

        print(f'scores for components {r} and {c}:', roc_aucs[r, c], accuracies[r, c])

print('-------- summary --------')
print('roc_aucs:')
print(roc_aucs)
print('accuracies:')
print(accuracies)

mask = np.triu(np.ones_like(accuracies, dtype=bool), 1)


# plot simple heatmap first
heat = plt.figure()
ax1 = heat.add_subplot(1, 1, 1) # heat.add_subplot(1, 2, 1)
# ax2 = heat.add_subplot(1, 2, 2)

sns.heatmap(
    accuracies, 
    mask=mask,
    ax=ax1, 
    xticklabels=range(1, num_components + 1), 
    yticklabels=range(1, num_components + 1), 
    cmap='magma',
    center=0.75, 
    vmin=0.5,
    vmax=1,
    cbar_kws={'label': 'prediction accuracy'}
)

# sns.heatmap(
#         roc_aucs, 
#         ax=ax2, 
#         xticklabels=range(1, num_components + 1), 
#         yticklabels=range(1, num_components + 1), 
#         cmap='magma',
#         center=0, 
#         vmin=0,
#         vmax=1,
#         cbar_kws={'label': 'roc auc'}
#     )

ax1.set_xlabel('Component')
ax1.set_ylabel('Component')

# ax2.set_xlabel('Component')
# ax2.set_ylabel('Component')
# ax2.set_title('ROC AUC')

heat.savefig(fig_path + '/triangle.png')