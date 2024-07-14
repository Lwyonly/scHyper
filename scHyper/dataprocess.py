import pandas as pd
import scanpy as sc
import numpy as np
import torch
import scipy
import matplotlib
import matplotlib.pyplot as plt
from itertools import product
import pkg_resources
import scHyper as ch

def trim_mean(arr, trim):

    sorted_arr = np.sort(arr)
    trim_count = int(len(arr) * trim)
    if trim_count == 0:
        trimmed_arr = sorted_arr
    else:
        trimmed_arr = sorted_arr[trim_count:-trim_count]
    if len(trimmed_arr) == 0:
        return 0
    else:
        return np.mean(trimmed_arr)

def meanExpression(adata, groupby, type="trimean", trim=0.25, return_df=False, normalization=False):

    from scipy.sparse import issparse
    from sklearn.preprocessing import normalize

    if scipy.sparse.issparse(adata.X):
        data = adata.X.todense()
        idx = adata.var_names
        col = adata.obs_names
        df = pd.DataFrame(data.T, index=idx, columns=col)
    else:
        data = adata.X.T
        idx = adata.var_names
        col = adata.obs_names
        df = pd.DataFrame(data, index=idx, columns=col)

    if type == "trimean":
        result = df.T.groupby(by=adata.obs[groupby]).apply(lambda group: group.apply(
        lambda x: trim_mean(x, trim), axis=0))
    else:
        result = df.T.groupby(by=adata.obs[groupby]).mean()

    if normalization:
        df_norm = normalize(result, copy=True, axis=1, norm="l2")
        result = pd.DataFrame(df_norm, columns=result.columns, index=result.index)

    if return_df is True:
        return result
    else:
        adata.uns["gene_call"] = result.T.to_dict()
        return adata

def process_ligands_receptors(adata, species_option, highly_variable=True):

    file_path = pkg_resources.resource_filename(
        __name__, f"database/{species_option}.csv")
    ligand_receptor_data = pd.read_csv(file_path,  header=0)
    ligand_receptor_data = ligand_receptor_data.set_index('interaction_name')
    ligand_receptor_data['interaction_name'] = ligand_receptor_data.index

    if highly_variable:
        def extract_genes(row):
            genes = row.split('_')
            return all(gene in highly_variable_genes for gene in genes)

        highly_variable_genes = adata.var[adata.var['highly_variable']].index
        ligand_receptor_data = ligand_receptor_data[ligand_receptor_data['Ligand_symbol'].apply(extract_genes)]
        ligand_receptor_data = ligand_receptor_data[ligand_receptor_data['Receptor_symbol'].apply(extract_genes)]
    else:
        ligand_receptor_data = ligand_receptor_data[ligand_receptor_data['Ligand_symbol'].isin(adata.var_names)]
        ligand_receptor_data = ligand_receptor_data[ligand_receptor_data['Receptor_symbol'].isin(adata.var_names)]

    gene_call = pd.DataFrame(adata.uns["gene_call"])
    unique_ligands = pd.Series(ligand_receptor_data['Ligand_symbol'].unique())
    unique_receptors = pd.Series(ligand_receptor_data['Receptor_symbol'].unique())
    gene_call_index = pd.Index(gene_call.index)

    split_genes_ligands = unique_ligands.str.split('_')
    exploded_genes_ligands = split_genes_ligands.explode()
    split_genes_receptors = unique_receptors.str.split('_')
    exploded_genes_receptors = split_genes_receptors.explode()

    common_ligand_genes = exploded_genes_ligands[exploded_genes_ligands.isin(gene_call_index)]
    common_receptor_genes = exploded_genes_receptors[exploded_genes_receptors.isin(gene_call_index)]
    ligands_gene_call = gene_call.loc[common_ligand_genes]
    receptors_gene_call = gene_call.loc[common_receptor_genes]
    adata.uns["ligands"] = ligands_gene_call.to_dict()
    adata.uns["receptors"] = receptors_gene_call.to_dict()

    return adata, ligand_receptor_data

def generate_tensor(adata, ligand_receptor_data):
    ligands_gene_call = pd.DataFrame(adata.uns["ligands"])
    receptors_gene_call = pd.DataFrame(adata.uns["receptors"])

    #ligands_data = ligands_gene_call.values
    #receptors_data = receptors_gene_call.values
    ligand_receptor_pairs = ligand_receptor_data.iloc[:,0:2].values

    tensor_shape = (len(ligands_gene_call.columns), len(receptors_gene_call.columns), len(ligand_receptor_data))
    tensor = np.zeros(tensor_shape)

    for i, (ligand, receptor) in enumerate(ligand_receptor_pairs):
        ligand_genes = ligand.split('_')
        ligand_nroot = np.power(np.prod(
            [ligands_gene_call.loc[gene,:].values if gene in ligands_gene_call.index else 1 for gene in
             ligand_genes], axis=0), 1 / len(ligand_genes))

        receptor_genes = receptor.split('_')
        receptor_nroot = np.power(np.prod(
            [receptors_gene_call.loc[gene,:].values if gene in receptors_gene_call.index else 1 for gene in
             receptor_genes], axis=0), 1 / len(receptor_genes))
        tensor[:, :, i] = np.outer(ligand_nroot, receptor_nroot)

    interaction_tensor = tensor.transpose(0, 2, 1)
    return interaction_tensor

def generate_triplets_weights_validlrindices(interaction_tensor):

    I, J, K = interaction_tensor.shape
    triplets = []
    weights = []
    validlrindices = []
    for i in range(I):
        for j in range(J):
            has_non_all_zero = 0
            for k in range(K):
                element = interaction_tensor[i, j, k]
                if element != 0:
                    triplets.append((i, j, k))
                    weights.append(element)
                    has_non_all_zero += 1
            if has_non_all_zero != 0:
                if j not in validlrindices:
                    validlrindices.append(j)
    triplets = np.array(triplets)
    weights = np.array(weights)
    validlrindices = np.array(validlrindices)

    return triplets, weights, validlrindices

def generate_validsenderindices_validreceiverindices(interaction_tensor):

    interaction_tensor = interaction_tensor.transpose(1, 0, 2)
    I, J, K = interaction_tensor.shape
    validsenderindices = []
    for i in range(I):
        for j in range(J):
            has_non_all_zero = 0
            for k in range(K):
                element = interaction_tensor[i, j, k]
                if element != 0:
                    has_non_all_zero += 1
            if has_non_all_zero != 0:
                if j not in validsenderindices:
                    validsenderindices.append(j)
    validsenderindices = np.array(validsenderindices)

    interaction_tensor = interaction_tensor.transpose(0, 2, 1)
    I, J, K = interaction_tensor.shape
    validreceiverindices = []
    for i in range(I):
        for j in range(J):
            has_non_all_zero = 0
            for k in range(K):
                element = interaction_tensor[i, j, k]
                if element != 0:
                    has_non_all_zero += 1
            if has_non_all_zero != 0:
                if j not in validreceiverindices:
                    validreceiverindices.append(j)
    validreceiverindices = np.array(validreceiverindices)

    return validsenderindices, validreceiverindices

def generate_validlrs_invalidlrs(validlrindices, ligand_receptor_data):

    df = pd.DataFrame({'j_indices': validlrindices})
    df = df.sort_values(by='j_indices')
    df = df.reset_index(drop=True)
    rows_to_extract = df['j_indices'].tolist()
    ligand_receptor_data.reset_index(drop=True, inplace=True)
    validlrs = ligand_receptor_data.iloc[rows_to_extract]
    invalidlrs = ligand_receptor_data.loc[~ligand_receptor_data.index.isin(rows_to_extract)]
    validlrs.reset_index(drop=True, inplace=True)
    invalidlrs.reset_index(drop=True, inplace=True)

    return validlrs, invalidlrs

def generate_validsenders_validreceivers(adata, validsenderindices, validreceiverindices):

    df = pd.DataFrame({'i_indices': validsenderindices})
    df = df.sort_values(by='i_indices')
    df = df.reset_index(drop=True)
    rows_to_extract = df['i_indices'].tolist()
    celltype = pd.DataFrame(adata.uns["gene_call"]).columns
    celltype = pd.DataFrame(celltype)
    validsenders = celltype.iloc[rows_to_extract]
    invalidsenders = celltype.iloc[~celltype.index.isin(rows_to_extract)]

    df_1 = pd.DataFrame({'k_indices': validsenderindices})
    df_1 = df_1.sort_values(by='k_indices')
    df_1 = df_1.reset_index(drop=True)
    rows_to_extract_1 = df_1['k_indices'].tolist()
    validreceivers = celltype.iloc[rows_to_extract_1]
    invalidreceivers = celltype.iloc[~celltype.index.isin(rows_to_extract_1)]

    return validsenders, invalidsenders, validreceivers, invalidreceivers

def update_weights(weights):
    max_weight = np.max(weights)
    min_weight = np.min(weights)
    scaling_factor = 5/ min_weight
    weights = weights * scaling_factor

    return weights

def update_triplets(triplets):

    for col in range(triplets.shape[1]):
        original_col = triplets[:, col]
        unique_nodes = np.unique(original_col)
        sorted_nodes = np.sort(unique_nodes)
        col_mapping = {node: new_node for new_node, node in enumerate(sorted_nodes)}
        remapped_col = np.vectorize(col_mapping.get)(original_col)
        triplets[:, col] = remapped_col

    return triplets

def generate_train_test(triplets, weights):

    random_indices = np.arange(len(triplets))
    np.random.shuffle(random_indices)
    triplets = triplets[random_indices]
    weights = weights[random_indices]
    train_triplets = np.array([]).reshape(0, 3)
    test_triplets = np.array([]).reshape(0, 3)
    train_weights = np.array([]).reshape(0,)
    test_weights = np.array([]).reshape(0,)
    unique_nodes = np.unique(triplets[:, 1])
    for node in unique_nodes:
        node_samples = triplets[triplets[:, 1] == node]
        num_samples = len(node_samples)
        if num_samples == 1:
            num_train_samples = 1
        else:
            num_train_samples = int(0.8 * num_samples)
        train_samples = node_samples[:num_train_samples]
        test_samples = node_samples[num_train_samples:]
        train_weights_samples = weights[np.where(triplets[:, 1] == node)][:num_train_samples]
        train_weights = np.hstack((train_weights, train_weights_samples))
        test_weights_samples = weights[np.where(triplets[:, 1] == node)][num_train_samples:]
        test_weights = np.hstack((test_weights, test_weights_samples))
        train_triplets = np.vstack((train_triplets, train_samples))
        test_triplets = np.vstack((test_triplets, test_samples))

    return train_triplets, test_triplets, train_weights, test_weights

def generate_nums_type(train_triplets, test_triplets):

    train_nums_type = [np.size(np.unique(train_triplets[:, 0])), np.size(np.unique(train_triplets[:, 1])), np.size(np.unique(train_triplets[:, 2]))]
    test_nums_type = [np.size(np.unique(test_triplets[:, 0])), np.size(np.unique(test_triplets[:, 1])), np.size(np.unique(test_triplets[:, 2]))]

    return train_nums_type, test_nums_type

def use_to_predict(triplets):
    unique_nodes_1 = np.unique(triplets[:, 0]).size
    unique_nodes_2 = np.unique(triplets[:, 1]).size
    unique_nodes_3 = np.unique(triplets[:, 2]).size
    use_to_predict = []
    for a in range(unique_nodes_1):
        for b in range(unique_nodes_2):
            for c in range(unique_nodes_3):
                use_to_triplets = (a,b,c)
                use_to_predict.append(use_to_triplets)
    use_to_predict = np.array(use_to_predict)
    return use_to_predict

def find_pattern_position(data, pattern):

    positions = []
    for subpattern in pattern:
        subpattern_array = np.array(subpattern)
        found_indices = np.where(np.all(data == subpattern_array, axis=1))[0]
        positions.append(found_indices)

    return positions


def genenrate_df_nn_candidates(validlrs, validsenders, validreceivers, triplets, use_to_predict):

    import pickle
    file_path = pkg_resources.resource_filename(
        __name__, f"model/predictions.pkl")
    with open(file_path, 'rb') as file:
        predictions = pickle.load(file)

    predictions = predictions.cpu()
    df_nn = pd.DataFrame(predictions.numpy(), columns=['Column1'])
    df_nn.rename(columns={'Column1': 'prob'}, inplace=True)

    pairs = [validsenders.iloc[:, 0].tolist(), validlrs['interaction_name'].tolist(), validreceivers.iloc[:, 0].tolist()]
    CT_LR_CT_pairs = [f'{pair1}_{pair2}_{pair3}' for pair1, pair2, pair3 in product(*pairs)]
    df_nn.index = CT_LR_CT_pairs

    matching_indices = find_pattern_position(use_to_predict, triplets)
    candidates = [CT_LR_CT_pairs[i] for indices_list in matching_indices for i in indices_list]

    return df_nn, candidates


def null_test(df_nn: pd.DataFrame, candidates, pval=0.05, plot=False):

    prob_test = df_nn[df_nn.index.isin(candidates)].copy()
    prob_null = df_nn[~df_nn.index.isin(candidates)]
    prob_test['p_val'] = prob_test['prob'].apply(
        lambda x: 1 - scipy.stats.percentileofscore(prob_null['prob'], x) / 100)
    df_enriched = prob_test[prob_test['p_val'] < pval].sort_values(by=['prob'], ascending=False)
    print(f'\nTotal enriched: {len(df_enriched)} / {len(df_nn)}')
    df_enriched['enriched_rank'] = np.arange(len(df_enriched)) + 1
    tensor_pval = torch.ones_like(torch.tensor(df_nn['prob'].values, dtype=torch.float32).view(-1, 1))
    tensor_pval[df_nn.index.isin(candidates), :] = torch.tensor(prob_test['p_val'].values, dtype=torch.float32).view(-1,1)

    if plot:
        cut = np.percentile(prob_null['prob'].values, pval)  # left tail
        plt.hist(prob_null['prob'], bins=1000, color='royalblue')
        for d in prob_test['prob']:
            if d < cut:
                c = 'red'
            else:
                c = 'gray'
            plt.axvline(d, ls=(0, (1, 1)), linewidth=0.5, alpha=0.8, c=c)
        plt.xlabel('prob')
        plt.show()
    del prob_test

    return df_enriched, tensor_pval

def generate_data_signaling(adata_h, file_path):

    import os
    ligands = pd.DataFrame(adata_h.uns["ligands"])
    receptors = pd.DataFrame(adata_h.uns["receptors"])
    ligands_receptors = pd.concat([ligands, receptors], axis=0)
    ligands_receptors = ligands_receptors.drop_duplicates()
    X = pd.DataFrame(adata_h.X)
    X.index = adata_h.obs.index
    X.columns = adata_h.var.index
    matching_columns = X.columns.intersection(ligands_receptors.index)
    data_signaling = X[matching_columns].T
    data_signaling.to_csv(os.path.join(file_path, 'data_signaling.csv'), index=True, header=True)

    return data_signaling

def generate_results(adata, df_nn, tensor_pval, validlrs, train_nums_type, file_path):

    import os
    interaction_df = pd.DataFrame(np.array(df_nn).reshape(-1, train_nums_type[0] * train_nums_type[1]))
    celltype = pd.DataFrame(adata.uns["gene_call"]).columns
    tensor_pval = pd.DataFrame(tensor_pval.reshape(-1, train_nums_type[0] * train_nums_type[1]))
    cell_type_counts = adata.obs['labels'].value_counts()

    interaction_df.to_csv(os.path.join(file_path, 'interaction_tensor.csv'), index=False, header=False)
    tensor_pval.to_csv(os.path.join(file_path, 'tensor_pval.csv'), index=False, header=False)
    adata.obs.to_csv(os.path.join(file_path, 'meta.csv'), index=True, header=True)
    validlrs.to_csv(os.path.join(file_path, 'validlrs.csv'), index=False, header=True)
    data_signaling = generate_data_signaling(adata, file_path=file_path)
