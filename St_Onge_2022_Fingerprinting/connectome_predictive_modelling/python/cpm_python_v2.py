"""
# Connectome Predictive Modelling - Python adaptation - V2.0

Script - CPM (LOOCV)
Author: Frederic St-Onge (adapted from Shen et al. 2017)

Date created: 2022-02-21
Date modified: 2022-11-03

Version 2.0

Future versions:
- Correlations in Matlab go much faster. We did some testing where we imported the original matlab function in Python. It works, but I am not sure if the gain in complexity is warranted for a little bit of time (really not that much time) gained.
- Parallelization with numba?
- Need to restore the other options from the original CPM (other types of feature selections, etc.)

"""

#CPM - Python translation
import os
import re
from argparse import ArgumentParser
import warnings


import numpy as np
import pandas as pd
import scipy
from sklearn import model_selection
from sklearn.linear_model import LinearRegression as lm
from sklearn.svm import SVR
from sklearn import metrics
from sklearn import preprocessing


def main():
    """ Main function launching the CPM.
    """
    print('---------------------------------------------')
    print('--Connectome Predictive Modelling in Python--')
    args = parsing()
    print('Parameters:')
    print(f'    Path to connectivity matrices (modality 1): {args.path_arrays_mod1}')
    if args.path_arrays_mod2:
        print(f'    Path to connectivity matrices (modality 2): {args.path_arrays_mod2}')
    print(f'    File to use for behavior predictions: {args.data_to_pred_path}')
    print('    ')
    print(f'    Type of model to train: {args.model}')
    print(f'    Threshold level for feature selection: {args.thresh_corr}')
    print(f'    Threshold level for shared edges in final mask: {args.thresh_retain}')
    print(f'    Path to output: {args.output}')
    print(f'    Name of the run: {args.name}')
    print(f'    ')
    print('---------------------------------------------')
    print(' ')

    print('Importing the arrays for modality 1...')
    behav_data, dict_id_arrays1 = dictionary_importer(args.path_arrays_mod1, args.data_to_pred_path)
    if args.path_arrays_mod2:
        print('Importing the arrays for modality 2...')
        dict_id_arrays2 = dictionary_importer(args.path_arrays_mod2)
    else:
        dict_id_arrays2 = None

    print('Merging the connectivity matrices together...')
    matrix_3d = matrix_progenitor(dict_id_arrays_mod1=dict_id_arrays1,
        dict_id_arrays_mod2=dict_id_arrays2, behav_data=behav_data)

    print('Launching CPMPY...')
    behav_train, behav_pred_train, behav_test, dict_rmse, dict_behav_pred = cpmpy(
        behav_data=behav_data, matrix_3d=matrix_3d, args=args)

    print('Computing CPM performance...')
    predict_meas_valid = cpm_prediction_measures(behav_train=behav_train,
        behav_pred_train=behav_pred_train, behav_test=behav_test, 
        dict_behav_pred=dict_behav_pred, args=args)

    print('Exporting CPM performance to file...')
    cpm_export(predict_meas_valid=predict_meas_valid, args=args)

def parsing():
    """ Take user arguments and return them.
    """

    parser = ArgumentParser(description="")

    parser.add_argument("-p1", "--path_arrays_mod1", help="Path to the connectivity matrices to use.")
    parser.add_argument("-p2", "--path_arrays_mod2", default=None, help="Path to connectivity matrices to the second modality to use.")
    parser.add_argument("-d", "--data_to_pred_path", help="Behavior data to predict with the connectivity.")
    parser.add_argument("-t", "--thresh_corr", default=0.01, type=float, help='In the CPM, threshold at which a p-value is considered as positive (at the edge-level).')
    parser.add_argument("-r", '--thresh_retain', default=0.95, type=float, help='Percentage of edges that should be retained in the final set.')
    parser.add_argument("-o", "--output", help="Path where the output should be stored once run.")
    parser.add_argument("-m", "--model", choices=['LM', 'SVR'], help='Which prediction model to use? If SVR, a grid search is performed.')
    parser.add_argument("-n", "--name", help="Suffix to add to the files to save.")

    args = parser.parse_args()

    return args

def dictionary_importer(path_arrays, data_to_pred_path=None):
    """ Imports the arrays and the IDs to analyse using a dictionary.

    Needs to run as many times as I have modalities.
    """

    #Import array names and ID data.
    folder_ls = os.listdir(path_arrays)
    behav_data = pd.read_csv(data_to_pred_path, index_col=0)

    #Empty dictionary to store the connectivity arrays
    dict_id_arrays = {}

    for file in folder_ls:
        #Extract ID from the files in the folder
        #TODO: Change below to an argument given by user for subject ID
        ids = re.search(r'[0-9]{6}', file).group()[0:6]

        #If the ID is in the target behavioral data, we import the corresponding array and store it
        if ids in behav_data.index.astype(str):
            array = np.loadtxt(f"{path_arrays}/{file}")

            dict_id_arrays[ids] = array

    if data_to_pred_path:
        return behav_data, dict_id_arrays
    else:
        return dict_id_arrays

def _3d_mat_baker(dicto):
    """ Hidden function. Takes a dictionary and turn the values in a 3d matrix.
    """
    list_arrays = []
    for arrays in dicto.values():
        list_arrays.append(arrays)

    matrix_3d = np.dstack(list_arrays)

    return matrix_3d

def matrix_progenitor(dict_id_arrays_mod1, dict_id_arrays_mod2, behav_data):
    """ Function taking a dictionary of ids and arrays and returning a 3D array instead.

    When we have 2 modalities, we concatenate them on the first axis (stick the cubes side-by-side)
    """

    matrix_3d_mod1 = _3d_mat_baker(dict_id_arrays_mod1)

    if dict_id_arrays_mod2:
        matrix_3d_mod2 = _3d_mat_baker(dict_id_arrays_mod2)
        matrix_3d = np.concatenate((matrix_3d_mod1, matrix_3d_mod2), axis=1)
    else:
        matrix_3d = matrix_3d_mod1.copy()

    #Run a couple user checks
    ### Check the number of matrices match the file inputed
    if matrix_3d.shape[2] == len(behav_data.index.values):
        print(f'    Total of {matrix_3d.shape[2]} participants')
    else:
        raise RuntimeError("    ERROR: The number of arrays doesn't match the desired"
        " number of participants")

    ### If second modality given, then should be same size as first
    if dict_id_arrays_mod2:
        if matrix_3d.shape[1] != (2 * matrix_3d.shape[0]):
            raise RuntimeError("    ERROR: The second modality"
            " should have the same number of nodes/columns as the first.")
    ### If only one modality, then should be symmetrical
    else:
        if matrix_3d.shape[0] != matrix_3d.shape[1]:
            raise RuntimeError("    ERROR: Currently, only "
            "symetric matrices are supported by cpmpy")

    return matrix_3d

def _cpm_prep_sizes(matrix_3d):
    """ Simple function computing size of arrays for the CPM.
    """

    nb_nodes_rows = matrix_3d.shape[0]
    nb_nodes_cols = matrix_3d.shape[1]

    return nb_nodes_rows, nb_nodes_cols

def _cpm_train_test_split(behav_data, matrix_3d):
    """ Wrapper function on top of sklearn's model selection.
    It just adds some data manipulation before and after to play
    around the 3d matrix.

    Returns the 2D split data
    """

    #Turn the FP values to a single flat numpy array
    behav_array = np.squeeze(behav_data.to_numpy())

    #Turn the 3D matrix to a 2D array where rows are edges and columns are IDs.
    flat_fc_array = matrix_3d\
        .reshape((matrix_3d.shape[0] * matrix_3d.shape[1]), matrix_3d.shape[2])\
        .T #Flip so the IDs are rows (samples)

    #Train test split
    mat_train, mat_test, behav_train, behav_test = model_selection\
        .train_test_split(flat_fc_array, behav_array, test_size=0.15, random_state=667)

    return mat_train, mat_test, behav_train, behav_test

def _cpm_edges_cube_reshape(mat, nb_nodes_rows, nb_nodes_cols):
    """ Short function reshaping the train and test matrices to cubes. Currently necessary to match
    the original matlab code, but can probably change later
    """

    #We use a Fortran order because we shifted the orientation of the dataset
    #  with a transpose earlier

    cube_mat = mat.reshape(nb_nodes_rows, nb_nodes_cols, mat.shape[0], order="F")

    return cube_mat

def _cpm_prep_empty(nb_nodes_rows, nb_nodes_cols, mat_train, mat_test):
    """ From sizes computed, returns empty arrays we will use to store information.
    """

    behav_pred_train = {}
    behav_pred_train["norm_pos"] = np.zeros((mat_train.shape[0], 1))
    behav_pred_train["norm_neg"] = np.zeros((mat_train.shape[0], 1))

    behav_pred_test = {}
    behav_pred_test["norm_pos"] = np.zeros((mat_test.shape[0], 1))
    behav_pred_test["norm_neg"] = np.zeros((mat_test.shape[0], 1))

    edges_cv_cube = {}
    edges_cv_cube['norm_pos'] = np.zeros((nb_nodes_rows, nb_nodes_cols,  mat_train.shape[0]))
    edges_cv_cube['norm_neg'] = np.zeros((nb_nodes_rows, nb_nodes_cols, mat_train.shape[0]))

    rmse_train_final = {}
    rmse_train_final['norm_pos'] = pd.DataFrame(columns=['c_value',
        'g_value', 'ep_value', 'rmse_pos', 'rmse_neg'])
    rmse_train_final['norm_neg'] = pd.DataFrame(columns=['c_value',
        'g_value', 'ep_value', 'rmse_pos', 'rmse_neg'])

    return behav_pred_train, behav_pred_test, edges_cv_cube, rmse_train_final

def _cpm_loocv_sample(mat_train, behav_train, leftout):
    """ Function iteratively leaving out one participant and returning the sample to 
    use for LOOCV
    """

    mat_train_loocv = np.delete(mat_train, leftout, axis=0)
    behav_train_loocv = np.delete(behav_train, leftout, axis=0)

    return mat_train_loocv, behav_train_loocv

def _cpm_feature_selection(mat_train_loocv, behav_train_loocv, nb_nodes_rows, nb_nodes_cols):
    """ Core of the CPM. Each edge is assessed for it's correlation with the fingerprinting.
    """

    r_val_arr = np.array([])
    p_val_arr = np.array([])

    s_behav_train_loocv = np.squeeze(behav_train_loocv)

    ########### Python-style correlation
    for col in range(mat_train_loocv.shape[1]):
        if (((col + 1) == 1) or ((col + 1) == mat_train_loocv.shape[1])):
            print(f'    ...Starting correlations {col + 1}/{mat_train_loocv.shape[1]}')

        edge_values = np.squeeze(mat_train_loocv[:,col])

        with warnings.catch_warnings():
            warnings.filterwarnings('error') #If warning, consider as error
            try:
                r_val, p_val = scipy.stats.pearsonr(x=edge_values, y=s_behav_train_loocv)
            except Warning: #If warning, set correlation to null
                r_val = 0
                p_val = 1

        r_val_arr = np.append(r_val_arr, r_val)
        p_val_arr = np.append(p_val_arr, p_val)

    r_val_mat = np.reshape(r_val_arr, (nb_nodes_rows, nb_nodes_cols))
    p_val_mat = np.reshape(p_val_arr, (nb_nodes_rows, nb_nodes_cols))

    return r_val_mat, p_val_mat

def _cpm_significant_edges(r_val_mat, p_val_mat, thresh_corr):
    """ Based on the correlation and p-value, we determine which edge is significant.
    """

    pos_mask = np.where(((r_val_mat > 0) & (p_val_mat < thresh_corr)), 1, 0)
    neg_mask = np.where(((r_val_mat < 0) & (p_val_mat < thresh_corr)), 1, 0)

    return pos_mask, neg_mask

def _cpm_sums(cube_mat, pos_mask, neg_mask, rand_pos_mask=None, rand_neg_mask=None):
    """ Function computing training sums from a cube.
    """

    sum_edges = {}
    sum_edges['norm_pos'] = np.zeros((cube_mat.shape[2], 1))
    sum_edges['norm_neg'] = np.zeros((cube_mat.shape[2], 1))

    try:
        if rand_pos_mask.any():
            sum_edges['rand_pos'] = np.zeros((cube_mat.shape[2], 1))
            sum_edges['rand_neg'] = np.zeros((cube_mat.shape[2], 1))
    except AttributeError:
        pass

    for sub in range(0, sum_edges['norm_pos'].shape[0]):
        sum_edges['norm_pos'][sub] = sum(sum(cube_mat[:,:,sub] * pos_mask)) / 2
        sum_edges['norm_neg'][sub] = sum(sum(cube_mat[:,:,sub] * neg_mask)) / 2
        #Check if array for random is empty.
        #Will throw error if yes, so we simply pass if we get error.
        try:
            if rand_pos_mask.any(): #Check if the array is empty
                sum_edges['rand_pos'][sub] = sum(sum(cube_mat[:,:,sub] * rand_pos_mask)) / 2
                sum_edges['rand_neg'][sub] = sum(sum(cube_mat[:,:,sub] * rand_neg_mask)) / 2
        except AttributeError:
            pass

    return sum_edges

def _cpm_test_loocv_sums(cube_mat, pos_mask, neg_mask, leftout):
    """ Function to compute the test sum during the loocv
    """

    test_mat = cube_mat[:,:,leftout]

    test_sum = {}
    test_sum["norm_pos"] = np.array([sum(sum(test_mat * pos_mask)) / 2])
    test_sum["norm_neg"] = np.array([sum(sum(test_mat * neg_mask)) / 2])

    return test_sum

def _svr_grid_search(train_sum_loocv, test_sum_loocv, behav_train_loocv, behav_test_loocv):
    """ Function performing a grid search on the SVR parameters. 
    """
    dict_rmse_loocv = {}

    param_data = pd.DataFrame(columns=['c_value', 'g_value',
        'ep_value', 'rmse_pos', 'rmse_neg'])

    for c_value in [1,10,100,1000,10000,100000]:
        for g_value in [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
            for ep_value in [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10]:

                behav_train_loocv_scaled, behav_test_loocv_scaled = _svr_scaler(behav_train_loocv,
                    behav_test_loocv)

                train_sum_loocv_scaled_pos, test_sum_loocv_scaled_pos = _svr_scaler(
                    train_sum_loocv['pos'],
                    test_sum_loocv['pos'])

                train_sum_loocv_scaled_neg, test_sum_loocv_scaled_neg = _svr_scaler(
                    train_sum_loocv['neg'],                    test_sum_loocv['neg'])

                rmse_test_pos_loocv = _svr_fit_predict(train_sum_loocv_scaled_pos,
                    test_sum_loocv_scaled_pos, behav_train_loocv_scaled,
                    behav_test_loocv_scaled, g_value=g_value, c_value=c_value, ep_value=ep_value, type_pred='loocv')

                rmse_test_neg_loocv = _svr_fit_predict(train_sum_loocv_scaled_neg,
                    test_sum_loocv_scaled_neg, behav_train_loocv_scaled,
                    behav_test_loocv_scaled, g_value=g_value, c_value=c_value, ep_value=ep_value, type_pred='loocv')

                params_final = {'c_value':[c_value], 'g_value':[g_value], 'ep_value':[ep_value],
                            'rmse_pos':[rmse_test_pos_loocv], 'rmse_neg':[rmse_test_neg_loocv]}

                params_final_df = pd.DataFrame.from_dict(params_final, orient="columns")

                param_data = pd.concat([param_data, params_final_df], axis=0, ignore_index=True)

    dict_rmse_loocv['norm_pos'] = param_data.loc[[param_data['rmse_pos'].idxmin()]]
    dict_rmse_loocv['norm_neg'] = param_data.loc[[param_data['rmse_neg'].idxmin()]]

    return dict_rmse_loocv

def _svr_scaler(train_sum, test_sum):
    """ Takes the imput to scale from train and test, scale them,
    seperate them, and return them.
    """

    sum_scale = np.concatenate((np.squeeze(train_sum), test_sum))
    sum_scale_res = sum_scale.reshape(-1,1)
    scaler_sum = preprocessing.StandardScaler()
    scaled_sum = scaler_sum.fit(sum_scale_res)\
        .transform(sum_scale_res)

    if test_sum.size == 1:
        train_sum_scaled = np.delete(scaled_sum, -1, axis=0)
        test_sum_scaled = scaled_sum[-1]
    else:
        train_sum_scaled = scaled_sum[0:train_sum.size]
        test_sum_scaled = scaled_sum[train_sum.size:]

    return train_sum_scaled, test_sum_scaled

def _svr_fit_predict(train_sum_scaled, test_sum_scaled, behav_train, behav_test,
    g_value, c_value, ep_value, type_pred):
    """ Fits and predicts and SVR model
    """
    svr_obj = SVR(kernel='rbf', gamma=g_value, C=c_value, epsilon=ep_value)
    svr_obj.fit(train_sum_scaled, behav_train.ravel())
    if type_pred == "final":
        behav_pred_test = svr_obj.predict(test_sum_scaled)
    else:
        behav_pred_test = svr_obj.predict(test_sum_scaled.reshape(1,-1))

    rmse_test = metrics.mean_squared_error(np.squeeze([behav_test]), np.squeeze([behav_pred_test]), squared=False)

    #When doing the prediction, if we do LOOCV we only return the RMSE.
    # When doing the prediction on the final sample, we return the predicted values as well
    if type_pred == "final":
        return rmse_test, behav_pred_test
    return rmse_test

def _svr_best_param_selection(rmse_train_final, type_m):
    """ Simple function grabbing the best parameters from the SVR grid search (lowest RMSE)
    """
    if type == "pos":
        min_rmse = rmse_train_final.loc[[rmse_train_final[f'rmse_{type_m}'].idxmin()]]
    else:
        min_rmse = rmse_train_final.loc[[rmse_train_final[f'rmse_{type_m}'].idxmin()]]

    c_value = min_rmse.loc[:,"c_value"].iloc[0]
    g_value = min_rmse.loc[:,"g_value"].iloc[0]
    ep_value = min_rmse.loc[:,"ep_value"].iloc[0]

    return c_value, g_value, ep_value

def _lm_fit_predict(train_sum, test_sum, behav_train, behav_test=None):
    """ Simple function fitting a predicting a linear model.
    When doing the final test prediction, we also want the RMSE.
    """

    fit_lm = lm().fit(train_sum, behav_train)

    try:
        if behav_test.any():
            behav_pred_test = fit_lm.predict(test_sum)
            rmse_test = metrics.mean_squared_error(np.squeeze([behav_test]), np.squeeze([behav_pred_test]), squared=False)

            return rmse_test, behav_pred_test
    except AttributeError:
        behav_pred_test = fit_lm.predict(test_sum.reshape(1,-1))
        return behav_pred_test

def _cpm_cross_validated_edges(edges_cv_cube, thresh_retain):
    """ Computes the final mask to use in the final train/test set
    """
    final_sum_edges = {}
    final_mask_edges = {}

    #Here, edges_cv_cube is a dictionary containing cubes of arrays, 
    # 1 for positive and 1 for negative edges
    for key, arrays in edges_cv_cube.items(): #For positive and negative edges...
        final_edges_cv = np.sum(arrays, axis=2)
        mask = np.where(final_edges_cv >= (thresh_retain*arrays.shape[2]), 1, 0)

        final_sum_edges[key] = final_edges_cv
        final_mask_edges[key] = mask

    return final_sum_edges, final_mask_edges

def _cpm_random_mask(pos_neg_masks):
    """ Here, we want to create two random masks of random edges for prediction.
    We create 2 because the number of edges might differ between the positive
    and negative edges. We want a random mask that has the same number of edges as the mask
    they are trying to mimic.

    Easiest way is to simply shuffle the original array
    """

    random_masks = {}

    for type_edges, mask in pos_neg_masks.items():
        flat_mask = mask.flatten() #Make the mask flat, because shuffle only works on 1 dimension
        np.random.seed(667) 
        np.random.shuffle(flat_mask) #Shuffle position of 0 and 1s (in-place)
        rand_mask = flat_mask.reshape(mask.shape) #Reshape the shuffled array to the original shape

        random_masks[type_edges] = rand_mask 

    return random_masks

def _cpm_export_edges(pos_neg_sums, args):
    """ Quick function exporting the matrices for the sum of participants
    using a given edge for predicted behavior.
    """
    for keys, arrays in pos_neg_sums.items():
        np.savetxt(f'{args.output}/sum_{keys}_edges_cv_{args.name}.csv', arrays, delimiter=',', 
            fmt='%f')

def cpmpy(behav_data, matrix_3d, args):
    """ Main CPM function. Under the hood it wraps multiple
    smaller functions that execute the different steps.
    """
    #First, we prepare the CPM
    ## Get information on sizes of matrices
    print("    Preparing CPM input...")
    nb_nodes_rows, nb_nodes_cols = _cpm_prep_sizes(matrix_3d=matrix_3d)

    ## Train-Test split of the data
    mat_train, mat_test, behav_train, behav_test = _cpm_train_test_split(behav_data=behav_data,
        matrix_3d=matrix_3d)

    ## Creating of empty arrays to store stuff during computation
    behav_pred_train, behav_pred_test, edges_cv_cube, rmse_train_final = _cpm_prep_empty(
        nb_nodes_rows=nb_nodes_rows, nb_nodes_cols=nb_nodes_cols,
        mat_train=mat_train, mat_test=mat_test)

    #Final step: for some operations, need a cube of the matrices. Reshape it here.
    cube_mat_train = _cpm_edges_cube_reshape(mat=mat_train,
        nb_nodes_rows=nb_nodes_rows, nb_nodes_cols=nb_nodes_cols)

    cube_mat_test = _cpm_edges_cube_reshape(mat=mat_test,
        nb_nodes_rows=nb_nodes_rows, nb_nodes_cols=nb_nodes_cols)

    # Now, we can launch the CPM
    print("    Launching cross-validation...")
    for leftout in range(0, mat_train.shape[0]):
        #if ((leftout + 1) == 1) or ((leftout + 1) % 20 == 0):
        print(f'    Leaving out participant {leftout + 1}/{mat_train.shape[0]}')

        mat_train_loocv, behav_train_loocv = _cpm_loocv_sample(mat_train=mat_train,
            behav_train=behav_train, leftout=leftout)

        r_val_mat, p_val_mat = _cpm_feature_selection(mat_train_loocv=mat_train_loocv,
            behav_train_loocv=behav_train_loocv, nb_nodes_rows=nb_nodes_rows,
            nb_nodes_cols=nb_nodes_cols)

        pos_mask, neg_mask =_cpm_significant_edges(r_val_mat=r_val_mat, p_val_mat=p_val_mat,
            thresh_corr=args.thresh_corr)

        edges_cv_cube["norm_pos"][:,:,leftout] = pos_mask
        edges_cv_cube["norm_neg"][:,:,leftout] = neg_mask

        #Might be a more elegant way of doing this, but I'm doing it this way
        #We basically need the matrices in a cube for easing some operations.

        cube_mat_loocv = np.delete(cube_mat_train, leftout, axis=2)

        train_sums_loocv = _cpm_sums(cube_mat_loocv, pos_mask=pos_mask, neg_mask=neg_mask)
        test_sums_loocv = _cpm_test_loocv_sums(cube_mat=cube_mat_train, pos_mask=pos_mask,
            neg_mask=neg_mask, leftout=leftout)

        behav_train_loocv = behav_train_loocv.copy()
        behav_test_loocv = np.array([behav_train[leftout]])

        if args.model == "LM":
            behav_pred_train['norm_pos'][leftout] = _lm_fit_predict(train_sums_loocv['norm_pos'],
                test_sums_loocv['norm_pos'], behav_train=behav_train_loocv)
            behav_pred_train['norm_neg'][leftout] = _lm_fit_predict(train_sums_loocv['norm_neg'],
                test_sums_loocv['norm_neg'], behav_train=behav_train_loocv)

        elif args.model == "SVR":
            svr_params, behav_pred_test[leftout] = _svr_grid_search(train_sums_loocv, test_sums_loocv, behav_train_loocv,
                behav_test_loocv)

            rmse_train_final['norm_pos'] = pd.concat([rmse_train_final['norm_pos'],
                svr_params['norm_pos']])
            rmse_train_final['norm_neg'] = pd.concat([rmse_train_final['norm_neg'],
                svr_params['norm_neg']])

    print("     Final validation in test set")
    #Final train/test measures for the final models
    final_sum_edges, final_mask_edges = _cpm_cross_validated_edges(
        edges_cv_cube=edges_cv_cube, thresh_retain=args.thresh_retain)

    #We write the final_sum_edges to file as we won't be needing it later and it's for future computations
    print(f'     ... Exporting cross-validated edges to file...')
    _cpm_export_edges(final_sum_edges, args=args)

    #Create a random mask that will be used to generate random sums for prediction.
    rando_mask = _cpm_random_mask(pos_neg_masks=final_mask_edges)

    #Compute all training/test sums
    train_sum_final = _cpm_sums(cube_mat=cube_mat_train, pos_mask=final_mask_edges['norm_pos'],
        neg_mask=final_mask_edges['norm_neg'], 
        rand_pos_mask=rando_mask['norm_pos'], rand_neg_mask=rando_mask['norm_neg'])

    test_sum_final = _cpm_sums(cube_mat=cube_mat_test, pos_mask=final_mask_edges['norm_pos'],
        neg_mask=final_mask_edges['norm_neg'], 
        rand_pos_mask=rando_mask['norm_pos'], rand_neg_mask=rando_mask['norm_neg'])

    #Final phase: Prediction metrics in final sample.
    if args.model == "SVR":
        dict_rmse = {} #For final RMSE
        dict_behav_pred = {} #For final prediction

        train_sum_final_scaled = {} #To store the scaled variables prior to prediction
        test_sum_final_scaled = {} 
        behav_scaled = {} 

        c_value_final = {} #Storing the hyper parameters for each type of run
        g_value_final = {}
        ep_value_final = {}

        #Looping over needed parameters, we scale the sums to use for SVR and we select the right
        #   svr parameters.
        for type_e in ['pos', 'neg']:
            for type_s in ['norm', 'rand']:
                train_sum_final_scaled[f'{type_s}_{type_e}'], 
                test_sum_final_scaled[f'{type_s}_{type_e}'] = _svr_scaler(
                    train_sum_final[f'{type_s}_{type_e}'], test_sum_final[f'{type_s}_{type_e}'])

                if type_s == "norm":
                    c_value_final[f'{type_s}_{type_e}'], g_value_final[f"{type_s}_{type_e}"],
                    ep_value_final[f'{type_s}_{type_e}'] = _svr_best_param_selection(
                        rmse_train_final=rmse_train_final[f'{type_s}_{type_e}'], type_m=f'{type_e}'
                    )

        #Last thing needed to scale is the behavior
        behav_scaled['train'], behav_scaled['test'] = _svr_scaler(
            behav_train, behav_test)

        #Next, we fit and predict our behavior using SVR
        for type_p in ['train', 'test']: #Type of prediction to do
            for type_e in ['pos', 'neg']: #Type of edge association (positive/negative)
                for type_s in ['norm', 'rand']: #Type of sum (normal or randomized?)
                    if type_p == "test":
                        test_sum_tmp = test_sum_final
                    elif type_p == "train":
                        test_sum_tmp = train_sum_final

                    dict_rmse[f'rmse_{type_s}_{type_p}_{type_e}'],
                    dict_behav_pred[f'behav_pred_{type_s}_{type_p}_{type_e}'] = _svr_fit_predict(
                        train_sum_scaled=train_sum_final_scaled[f'{type_s}_{type_e}'],
                        test_sum_scaled=test_sum_final_scaled[f'{type_s}_{type_e}'],
                        behav_train=behav_scaled[f'train'], behav_test=behav_scaled[f'{type_p}'],
                        c_value=c_value_final[f'{type_s}_{type_e}'], 
                        g_value=g_value_final[f"{type_s}_{type_e}"],
                        ep_value=ep_value_final[f'{type_s}_{type_e}'], type_pred='final')

    if args.model == "LM":
        dict_rmse = {}
        dict_behav_pred = {}

        #We iterate over the different metrics we need to extract, and store them in dictionaries
        for type_p in ['test', 'train']:
            for type_e in ['pos', 'neg']:
                for type_s in ['rand', 'norm']:
                    #To avoid having multiple calls to _lm_fit_predict, we use aliases for the
                    # variables to pass. Specifically, we want the prediction in both train and
                    # test set, so we just assign them to different aliases before running pred.
                    if type_p == "test":
                        test_sum_tmp = test_sum_final
                        behav_test_tmp = behav_test
                    elif type_p == "train":
                        test_sum_tmp = train_sum_final
                        behav_test_tmp = behav_train

                    dict_rmse[f'rmse_{type_s}_{type_p}_{type_e}'], dict_behav_pred[f'behav_pred_{type_s}_{type_p}_{type_e}'] = _lm_fit_predict(
                        train_sum=train_sum_final[f'{type_s}_{type_e}'],
                        test_sum=test_sum_tmp[f'{type_s}_{type_e}'],
                        behav_train=behav_train, behav_test=behav_test_tmp)

    #The CPM outputs many things:
        # 1) The actual behavior observed for the training and testing sample
        # 2) The predicted values during the cross-validation
        # 3) The actual behavior observed in the leftout test set
        # 4) The RMSE of the prediction in the leftout test set
        # 5) The arrays of predicted behavior for participants in the leftout test set

    return behav_train, behav_pred_train, behav_test, dict_rmse, dict_behav_pred

def cpm_prediction_measures(behav_train, behav_pred_train, behav_test, dict_behav_pred, args):
    """ Here, the idea is to compute the different prediction metrics we are interested in,
    and store that information for output.

    Metrics we output are:
    - For pos/neg edges
        - Cross-validation:
            - Correlation between original behavior and predicted behavior during LOOCV
            - RMSE between original behavior and predicted behavior during LOOCV
        - Leftout test set:
            - For normal and randomized FC sums
                - Correlation between test set behavior and predicted behavior (final)
                - RMSE between test set behavior and predicted behavior

    In total, we should have 12 different metrics
    """
    predict_meas_valid = pd.DataFrame(columns=['model', 'type_edge', 'type_sum', 'type_samp', 
        'type_predict', 'type_meas', 'value'])

    #First, compute for the cross-validation
    for type_e in ['pos', 'neg']:
        for type_s in ['norm', 'rand']:
            if type_s == 'norm':
                corr_val, p_val = scipy.stats.pearsonr(np.squeeze(behav_train),
                    np.squeeze(behav_pred_train[f'{type_s}_{type_e}']))
                rmse_val = metrics.mean_squared_error(y_true=np.squeeze(behav_train),
                    y_pred=np.squeeze(behav_pred_train[f'{type_s}_{type_e}']),
                        squared=False)

                tmp_df = pd.DataFrame(data={'model':[args.model, args.model, args.model], 
                'type_edge':[type_e, type_e, type_e],
                'type_sum':[type_s, type_s, type_s], 
                'type_samp':['cross_validation', 'cross_validation', 'cross_validation'],
                'type_predict':['norm', 'norm', 'norm'],
                'type_meas':['corr', 'pval', 'rmse'], 
                'value':[corr_val, p_val, rmse_val]})

                predict_meas_valid = pd.concat([predict_meas_valid, tmp_df], axis=0, 
                    ignore_index=True)

    #Next, compute for the test set
    for type_e in ['pos', 'neg']:
        for type_s in ['norm', 'rand']:
            for type_p in ['train', 'test']:
                if type_p == 'train':
                    corr_val, p_val = scipy.stats.pearsonr(np.squeeze(behav_train),
                        np.squeeze(dict_behav_pred[f'behav_pred_{type_s}_{type_p}_{type_e}']))
                    rmse_val = metrics.mean_squared_error(y_true=np.squeeze(behav_train),
                        y_pred=np.squeeze(dict_behav_pred[f'behav_pred_{type_s}_{type_p}_{type_e}']),
                            squared=False)
                elif type_p == "test":
                    corr_val, p_val = scipy.stats.pearsonr(np.squeeze(behav_test),
                        np.squeeze(dict_behav_pred[f'behav_pred_{type_s}_{type_p}_{type_e}']))
                    rmse_val = metrics.mean_squared_error(y_true=np.squeeze(behav_test),
                        y_pred=np.squeeze(dict_behav_pred[f'behav_pred_{type_s}_{type_p}_{type_e}']),
                            squared=False)

                tmp_df = pd.DataFrame(data={
                    'model':[args.model, args.model, args.model],
                    'type_edge':[type_e, type_e, type_e],
                    'type_sum':[type_s, type_s, type_s],
                    'type_samp':['final', 'final', 'final'], 
                    'type_predict':[type_p, type_p, type_p],
                    'type_meas':['corr', 'pval', 'rmse'],
                    'value':[corr_val, p_val, rmse_val]})

                predict_meas_valid = pd.concat([predict_meas_valid, tmp_df], axis=0, 
                    ignore_index=True)

    return predict_meas_valid

def cpm_export(predict_meas_valid, args):
    """ Simple function exporting the predicted measures
    """
    predict_meas_valid.to_csv(f"{args.output}/model_validation_{args.model}_{args.name}.csv")

if __name__ == "__main__":
	main()