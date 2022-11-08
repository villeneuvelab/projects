""" ICC analyses

First, threshold the entire matrix at 95th percentile. 

So now we have 4 matrices: the original ones, and the thresholded ones.

Results to compute:
- Overlap between windows (ICC or Jaccard Coefficient) at a whole brain level
- Proportion of edges by network + test within window of proportion

First step, split the matrix in networks by using a dictionary.


"""

import os
import re
import numpy as np
import pandas as pd
from sklearn.metrics import jaccard_score

#First step, import all the ICC data for a given window

def import_data_clean(path_to_data):
    """ Imports all data from a single window.
    """

    dict_edges_sub = {}
    dict_edges_task = {}

    #We iterate over the files
    for type_icc in ['sub', 'task']: #For both types of ICCs
        for file in os.listdir(path_to_data): #We iterate over the input path
            if type_icc in file:
                name = re.search(r'_w[0-9]{2}', file).group()[1:4]
                print(name)
                data = np.genfromtxt(f"{path_to_data}/{file}", delimiter=',')
                if type_icc == "sub":
                    dict_edges_sub[f'{name}'] = data
                else:
                    dict_edges_task[f'{name}'] = data

    return dict_edges_sub, dict_edges_task

def rois():
    """ Output dictionary of regions under consideration
    """
    dict_net_schaefer_within = {
        #Based on the dataset.labels from the Schaefer Atlas in Nilearn.
        "visual": list(range(0, 31)) + list(range(200, 230)),
        "somatomotor": list(range(31, 68)) + list(range(230, 270)),
        "dorsal_attention": list(range(68, 91)) + list(range(270, 293)),
        "salience_ventral_attention": list(range(91, 113)) + list(range(293, 318)),
        "limbic": list(range(113, 126)) + list(range(318, 331)),
        "frontoparietal": list(range(126, 148)) + list(range(331, 361)),
        "default_mode": list(range(148, 200)) + list(range(361, 400)),
        "whole_brain":list(range(0, 400))
        }

    #Dictionary for the Schaefer Atlas (Between- network) 
    dict_net_schaefer_between_all = {
        #Overall between-network connectivity (we don't do the average between each network)
        "visual": list(range(31, 200)) + list(range(231, 400)),
        "somatomotor": list(range(68, 200)) + list(range(270, 400)),
        "dorsal_attention": list(range(91, 200)) + list(range(293, 400)),
        "salience_ventral_attention": list(range(113, 200)) + list(range(318, 400)),
        "limbic": list(range(126, 200)) + list(range(331, 400)),
        "frontoparietal": list(range(148, 200)) + list(range(361, 400)),
        "default_mode": list(range(0, 148)) + list(range(200, 361))
        }

    return dict_net_schaefer_within, dict_net_schaefer_between_all

def bin_icc_vals(dict_edges):
    """ Dictionary of edges for which we need to return binary info.
    """
    dict_bin = {}
    dict_val_bin = {}

    for window, arrays in dict_edges.items():

        perc_val = np.nanpercentile(arrays, q=95, axis=None)

        print(perc_val)
        dict_bin[f'{window}'] = np.where(arrays >= perc_val, 1, 0)
        dict_val_bin[f'{window}'] = np.where(arrays >= perc_val, arrays, np.NaN)

    return dict_bin, dict_val_bin

def _mean_std_icc(submatrix):
    """ Takes a submatrix and computes the average
    """

    #If we have a symmetric matrix, we want average of just the upper, non-diagonal elements
    # for the metrics
    if submatrix.shape[0] == submatrix.shape[1]:
        triu_extract = submatrix[np.triu_indices_from(submatrix, k=1)]
        mean_val = np.nanmean(triu_extract)
        std_val = np.nanstd(triu_extract)

    else:
        mean_val = np.nanmean(submatrix)
        std_val = np.nanstd(submatrix)

    return mean_val, std_val

def _sum_binary(submatrix):
    """ Computes the sum of significant edges by network
    """

    if submatrix.shape[0] == submatrix.shape[1]:
        triu_extract = submatrix[np.triu_indices_from(submatrix,k=1)]
        sum_mat = np.nansum(triu_extract)
    else:
        sum_mat = np.nansum(submatrix)

    return sum_mat

def _prop_net(submatrix):
    """ Computes the proportion of binary edges by network.
    """

    if submatrix.shape[0] == submatrix.shape[1]:
        size_net = ((submatrix.shape[0] * submatrix.shape[1]) - len(submatrix)) / 2
        triu_extract = submatrix[np.triu_indices_from(submatrix, k=1)]
        prop_net = (np.sum(triu_extract) / size_net)
    else:
        size_net = submatrix.shape[0] * submatrix.shape[1]
        prop_net = (np.sum(submatrix) / size_net)

    return prop_net, size_net

def _icc_metrics_computer(pop, keys, type_meas, submatrix, type_mat):
    """ Determines what dataframe is output, depending on what is inputed.
    """

    if type_mat == "bin":
        sum_mat = _sum_binary(submatrix=submatrix)
        prop_mat, size_mat = _prop_net(submatrix=submatrix)

        tmp_df = pd.DataFrame(data={"window":[f"{pop}", f'{pop}',  f'{pop}'],
                    'network':[f'{keys}', f'{keys}', f'{keys}'],
                    'type':[f'{type_meas}', f'{type_meas}', f'{type_meas}'],
                    'metric':['sum_bin', 'prop_bin', 'size_network'],
                    'value':[sum_mat, prop_mat, size_mat]})

    elif type_mat == "icc":
        mean_val, std_val = _mean_std_icc(submatrix=submatrix)
        tmp_df = pd.DataFrame(data={"window":[f"{pop}", f'{pop}'],
                    'network':[f'{keys}', f'{keys}'],
                    'type':[f'{type_meas}', f'{type_meas}'],
                    'metric':['mean_icc', 'std_icc'],
                    'value':[mean_val, std_val]})

    return tmp_df

def icc_val_network_avg_bn(dict_nodes_within, dict_nodes_between, icc_dict, type_mat):
    """ Function returning a nested dictionary where we output the pure
    ICC values in each network.

    This function returns values between-network connections related to
    1 network (doesn't do each network-network edges)
    """
    icc_metrics = pd.DataFrame(columns=['window', 'network', 'type', 'metric', "value"])

    for pop, mat in icc_dict.items(): #For each window, we get the ICC matrix
        for keys1 in dict_nodes_within.keys(): #Iterate overt the networks
            #First, within-network edges
            for keys2 in dict_nodes_within.keys(): #Iterate again
                if keys1 == keys2: #Within network, so we only want the diagonal elements.
                    submatrix = mat[dict_nodes_within[f'{keys1}']][:, dict_nodes_within[f'{keys2}']]
                    type_meas = "within"

                    #Create a temporary dataframe with the right metrics computed.
                    tmp_df = _icc_metrics_computer(pop=pop, keys=keys1, type_meas=type_meas,
                        submatrix=submatrix, type_mat=type_mat)

                    #Store the metrics in the main dataframe
                    icc_metrics = pd.concat([icc_metrics, tmp_df],
                        ignore_index=True)

            #Second, between-network edges. We don't do each between network separately.
            for keys2 in dict_nodes_between.keys():
                if keys1 == keys2:
                    submatrix = mat[dict_nodes_within[f'{keys1}']][:,dict_nodes_between[f'{keys2}']]
                    type_meas = "between"

                    tmp_df = _icc_metrics_computer(pop=pop, keys=keys1, type_meas=type_meas,
                        submatrix=submatrix, type_mat=type_mat)

                    icc_metrics = pd.concat([icc_metrics, tmp_df],
                        ignore_index=True)

    return icc_metrics

def _reorder_dict(icc_bin_dict):
    """ Reorders the windows, based on their length, to make proper comparisons.
    """
    if len(icc_bin_dict) == 10:
        dict_order = ['w01', 'w02', 'w03', 'w04', 'w05', 'w06', 'w07', 'w08', 'w09', 'w10']
    elif len(icc_bin_dict) == 15:
        dict_order = ['w01', 'w02', 'w03', 'w04', 'w05', 'w06', 'w07', 'w08', 'w09', 'w10', 'w11', 'w12', 'w13', 'w14', 'w15']
    elif len(icc_bin_dict) == 13:
        dict_order = ['w01', 'w02', 'w03', 'w04', 'w05', 'w06', 'w07', 'w08', 'w09', 'w10', 'w11', 'w12', 'w13']
    elif len(icc_bin_dict) == 8:
        dict_order = ['w01', 'w02', 'w03', 'w04', 'w05', 'w06', 'w07', 'w08']

    reordered_dict = {k: icc_bin_dict[k] for k in dict_order}
    
    return reordered_dict

def icc_window_overlap_bin(dict_nodes_within, dict_nodes_between, icc_bin_dict):
    """ Goal is to compute and store the overlap between two windows in terms of 95th highest
    ICC values.
    """
    icc_overlap = pd.DataFrame(columns=['window_pair', 'network', 'type', 'value', 'size_net'])

    reordered_dict = _reorder_dict(icc_bin_dict)

    for net, values in dict_nodes_within.items():
        for num_1, (keys1, arrays1) in enumerate(reordered_dict.items()):
            submatrix1 = arrays1[values][:,values]
            fin_mat1 = submatrix1[np.triu_indices_from(submatrix1, k=1)]
            size_net = size_net = ((submatrix1.shape[0] * submatrix1.shape[1]) - len(submatrix1)) / 2

            for num_2, (keys2, arrays2) in enumerate(reordered_dict.items()):
                submatrix2 = arrays2[values][:,values]
                fin_mat2 = submatrix2[np.triu_indices_from(submatrix2, k=1)]

                if num_1 + 1 == num_2:
                    key_pair = keys1 + "_" + keys2
                    overlap = jaccard_score(np.ravel(fin_mat1), np.ravel(fin_mat2))

                    tmp_df = pd.DataFrame(data={"window_pair":[key_pair], "network":[net],
                    "type":["within"], "value":[overlap], 'size_net':[size_net]})

                    icc_overlap = pd.concat([icc_overlap, tmp_df])

    for net, values in dict_nodes_between.items():
        for num_1, (keys1, arrays1) in enumerate(reordered_dict.items()):
            fin_mat1 = arrays1[dict_nodes_within[net]][:,values]
            size_net = fin_mat1.shape[0] * fin_mat1.shape[1]

            for num_2, (keys2, arrays2) in enumerate(reordered_dict.items()):
                fin_mat2 = arrays2[dict_nodes_within[net]][:,values]

                if num_1 + 1 == num_2:
                    key_pair = keys1 + "_" + keys2
                    overlap = jaccard_score(np.ravel(fin_mat1), np.ravel(fin_mat2))

                    tmp_df = pd.DataFrame(data={"window_pair":[key_pair], "network":[net],
                    "type":["between"], "value":[overlap], 'size_net':[size_net]})

                    icc_overlap = pd.concat([icc_overlap, tmp_df])

    return icc_overlap.reset_index().set_index('window_pair').drop('index', axis=1)

def icc_sub_task_overlap(dict_nodes_within, dict_nodes_between, icc_bin_dict_sub,
    icc_bin_dict_task):
    """ Quick overlap computation between subject and task ICC.
    """
    icc_overlap_df = pd.DataFrame(columns=['window_pair', 'network', 'type',
        'val_overlap_sub_task', 'size_net'])

    for net, values in dict_nodes_within.items():
        for keys1, arrays1 in icc_bin_dict_sub.items():
            submatrix1 = arrays1[values][:,values]
            fin_mat1 = submatrix1[np.triu_indices_from(submatrix1, k=1)]
            size_net = ((submatrix1.shape[0] * submatrix1.shape[1]) - len(submatrix1)) / 2

            for keys2, arrays2 in icc_bin_dict_task.items():
                submatrix2 = arrays2[values][:,values]
                fin_mat2 = submatrix2[np.triu_indices_from(submatrix2, k=1)]
                if keys1 == keys2:
                    key_pair = keys1 + "_" + keys2
                    overlap = jaccard_score(np.ravel(fin_mat1), np.ravel(fin_mat2))

                    tmp_df = pd.DataFrame(data={'window_pair':[key_pair], 'network':[net],
                        'type':['within'], 'val_overlap_sub_task':[overlap], 'size_net':[size_net]})

                    icc_overlap_df = pd.concat([icc_overlap_df, tmp_df])

    for net, values in dict_nodes_between.items():
        for keys1, arrays1 in icc_bin_dict_sub.items():
            fin_mat1 = arrays1[dict_nodes_within[net]][:,values]
            size_net = fin_mat1.shape[0] * fin_mat1.shape[1]

            for keys2, arrays2 in icc_bin_dict_task.items():
                fin_mat2 = arrays2[dict_nodes_within[net]][:,values]
                if keys1 == keys2:
                    key_pair = keys1 + "_" + keys2
                    overlap = jaccard_score(np.ravel(fin_mat1), np.ravel(fin_mat2))

                    tmp_df = pd.DataFrame(data={'window_pair':[key_pair], 'network':[net],
                        'type':['between'], 'val_overlap_sub_task':[overlap], 'size_net':[size_net]})

                    icc_overlap_df = pd.concat([icc_overlap_df, tmp_df])

    return icc_overlap_df.reset_index().set_index('window_pair').drop('index', axis=1)

