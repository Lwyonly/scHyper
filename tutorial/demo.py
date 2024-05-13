import pandas as pd
import scanpy as sc
import numpy as np
import torch
import scipy
import matplotlib
import matplotlib.pyplot as plt
import os
from itertools import product
from scHyper import dataprocess


# Read gene expression matrix and cell type labels
x = pd.read_csv("to/path/tutorial/data/count.csv", index_col=0)
meta = pd.read_csv("to/path/tutorial/data/meta.csv", index_col=0)

# We add the labels to adata.obs, and please normalize and logarithmize the unprocessed data here.
adata_h = sc.AnnData(X=x.T.values)
adata_h.obs = pd.DataFrame(meta["labels"])
adata_h.var_names = x.index
adata_h.obs_names = x.columns
adata_h.obs_names_make_unique()
adata_h.var_names_make_unique()
sc.pp.filter_cells(adata_h, min_genes=10)
sc.pp.filter_genes(adata_h, min_cells=3)
sc.pp.normalize_total(adata_h, target_sum=1e4)
sc.pp.log1p(adata_h)

# Calculate the average expression of genes in different cell types
# Please select the truncated mean or average based on your requirements
adata_h = dataprocess.meanExpression(adata_h, type="mean", groupby="labels")
# Pair of ligand-receptor interactions, expression of ligands and receptors appear in the gene expression matrix
# Please select whether to use high-variation genes based on your needs
adata_h, ligand_receptor_data = dataprocess.process_ligands_receptors(adata_h, "human", highly_variable=False)
# Construct intercellular communication tensor
interaction_tensor = dataprocess.generate_tensor(adata_h, ligand_receptor_data)

# Get the triples of hypergraph and weights.
triplets, weights, validlrindices = dataprocess.generate_triplets_weights_validlrindices(interaction_tensor)
# Obtain effective L-R pairs and ineffective L-R pairs.
validlrs, invalidlrs = dataprocess.generate_validlrs_invalidlrs(validlrindices, ligand_receptor_data)
# Obtain effective celltypes and ineffective celltypes.
validsenderindices, validreceiverindices = dataprocess.generate_validsenderindices_validreceiverindices(interaction_tensor)
validsenders, invalidsenders, validreceivers, invalidreceivers = dataprocess.generate_validsenders_validreceivers(adata_h, validsenderindices, validreceiverindices)

# Update weights and triples
triplets = dataprocess.update_triplets(triplets)
weights = dataprocess.update_weights(weights)

# Generate the training set and test set
train_triplets, test_triplets, train_weights, test_weights = dataprocess.generate_train_test(triplets, weights)
train_nums_type, test_nums_type = dataprocess.generate_nums_type(train_triplets, test_triplets)

# Save the training and test data sets, Please specify the save_path
save_path='to/path/data/demo'
np.savez(os.path.join(save_path, 'train_data.npz'), nums_type=train_nums_type, train_data=train_triplets, train_weight=train_weights)
np.savez(os.path.join(save_path, 'test_data.npz'), nums_type=train_nums_type, test_data=test_triplets, test_weight=test_weights)

# Create and save an array for prediction
use_to_predict = dataprocess.use_to_predict(triplets)
np.savez(os.path.join(save_path, 'use_to_predict.npz'), use_to_predict=use_to_predict)

# The next step is after training the model
# We used nonparametric tests to identify significant intercellular communications
df_nn, candidates = dataprocess.genenrate_df_nn_candidates(validlrs, validsenders, validreceivers, triplets, use_to_predict)
df_enriched, tensor_pval = dataprocess.null_test(df_nn, candidates, pval=0.05, plot=False)

# Visualization preparation
file_path='to/path/results'
dataprocess.generate_results(adata_h, df_nn, tensor_pval, validlrs, train_nums_type, file_path)





