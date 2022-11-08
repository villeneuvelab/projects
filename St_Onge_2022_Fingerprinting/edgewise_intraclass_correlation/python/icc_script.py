""" ICC subject and task
"""

from argparse import ArgumentParser
import warnings
import re
import os
import numpy as np
import pandas as pd
import pingouin as pg

def main():
    """ Launching ICC analysis
    """
    print('---------------------------------------------')
    print('-------Edge-wise Intraclass Correlation------')
    args = parsing()
    print("Parameters:")
    print(f"    Path to the connectivity matrices (Mod 1): {args.path_arrays_1}")
    print(f"    Path to the connectivity matrices (Mod 2): {args.path_arrays_2}")
    print(f'    File containing the IDs: {args.path_ids}')
    print(f'    Path to the output: {args.output}')
    print(f'    Name of the run: {args.name}')
    print('---------------------------------------------')
    print(" ")
    print("Importing the arrays in a dictionary...")
    dict_arrays_1, dict_arrays_2 = dictionary_importer(args.path_arrays_1, path_arrays_2=args.path_arrays_2, path_ids=args.path_ids)
    print("Amalgamating the arrays...")
    concatenated_df = amalgamate(dict_arrays_1=dict_arrays_1, dict_arrays_2=dict_arrays_2)
    print("Computing Subject and Task Edge-wise ICC...")
    icc_sub_array, icc_task_array = edge_wise_icc(concatenated_df=concatenated_df)
    print("Reshaping ICC arrays to as a matrix...")
    res_icc_sub_array, res_icc_task_array = reshaping(icc_sub_array=icc_sub_array,
    icc_task_array=icc_task_array)
    print("Exporting resulting ICC matrices...")
    export(res_icc_sub_array=res_icc_sub_array,
    res_icc_task_array=res_icc_task_array, output=args.output, name=args.name)

def parsing():
    """ Parses the arguments given to the script
    """

    parser = ArgumentParser(description="")
    parser.add_argument("-p1", "--path_arrays_1", help='Path to the connectivity matrices of modality 1.')
    parser.add_argument("-p2", "--path_arrays_2", help='Path to the connectivity matrices of modality 2')
    parser.add_argument("-i", "--path_ids", help='Path to the spreadsheet with the IDs')
    parser.add_argument('-o', '--output', help='Path to where the paths should be output.')
    parser.add_argument('-n', '--name', help='Suffix to add to the files to save')

    args = parser.parse_args()

    return args

def dictionary_importer(path_arrays_1, path_arrays_2, path_ids):
    """ Function importing the arrays from both modalities to a single
    dictionary. 
    """

    #To store the connectivity matrices of both modalities
    dict_arrays_1 = {}
    dict_arrays_2 = {}

    #For each path (each modality)
    for i, paths in enumerate([path_arrays_1, path_arrays_2]):
        folder = os.listdir(paths) #List files
        id_df = pd.read_csv(path_ids) #Get IDs we want analyzed
        id_list = id_df.iloc[:,0].astype(str).to_list()

        for matrix in folder:
            #Use regex to extract the participant ID from the matrix file name
            id_part = re.search(r'[0-9]{6}', matrix).group()[0:6]

            #If the extracted ID is part of our list of subject of interest
            if id_part in id_list:
                array = np.loadtxt(f"{paths}/{matrix}") #Load the matrix
                #Store the matrix in the right dictionary
                if i == 0:
                    dict_arrays_1[id_part] = array
                elif i == 1:
                    dict_arrays_2[id_part] = array

    return dict_arrays_1, dict_arrays_2

def amalgamate(dict_arrays_1, dict_arrays_2):
    """ Transforms and congregates individual-level connectivity
    matrices to a single long-format dataframe. Each connectivity
    matrix is first flattened and used to create a long-format
    dataset for a single individual. Then, we concatenate the
    individual-level datasets in a single one. 
    """

    #To store individual-level data
    list_df = []

    #For each dictionary...
    for i, dicts in enumerate([dict_arrays_1, dict_arrays_2]):
        #We store the modality in use depending on which dict we are working
        if i == 0:
            mod = "rest"
        elif i == 1:
            mod = "task"

        #For each dictionary, we iterate the items (pairs of IDs and matrices)
        for keys, arrays in dicts.items():

            #Make the connectivity matrix flat
            flat_array = arrays.flatten()

            #Create new dataframe with 2 scalar (ID and modality) and the connectivity values
            tmp_df = pd.DataFrame(data={"id":keys,
            "mod":mod, "conn_value":flat_array})

            #To do the "edge-wise" part, we need to execute the next
            # operations on each edge individually. Here we just give a
            # sequential number to each row (i.e., each edge).
            tmp_df['edge'] = np.arange(len(tmp_df))

            list_df.append(tmp_df)
    
    #Merge all subjects together. Gives a 32 million rows dataset.
    concatenated_df = pd.concat(list_df, axis=0, ignore_index=True)

    return concatenated_df

def edge_wise_icc(concatenated_df):
    """ Computes the two ICC types for each edge individually.

    Note: in preliminary testing, when edges have constant values (e.g., 
    when all cells are 1 or 0, it throws a warning. Normally, it should be
    ignored, but it causes the system to crash, just like in the CPM).
    So now, edges that would cause a warning are simply replaced by missing.
    """

    unique_edges = concatenated_df['edge'].unique() #To iterate over edges (edge-wise)
    icc_sub_array = np.empty(unique_edges.shape) #To store the iccs
    icc_task_array = np.empty(unique_edges.shape)

    #For each individual edge...
    for i, edge in enumerate(unique_edges):
        #Restrict main dataset to the edge of interest
        restrict_df = concatenated_df[concatenated_df['edge'] == edge]

        if ((i + 1) == 1) or ((i + 1) % 10000 == 0):
            print(f"    Processing edge {i + 1}/{len(unique_edges)}")

        #Compute the ICCs
        with warnings.catch_warnings():
            warnings.filterwarnings('error') #If warning, consider as error
            try:
                #Subject-level ICC (Judge are IDs, targets are modalities)
                #Output is a dataframe.
                icc_sub_df = pg.intraclass_corr(data=restrict_df,
                targets="mod", raters='id', ratings="conn_value")\
                    .set_index('Type') #Setting index to type so can grab with loc
                #Extract just what we need from dataframe and save in array
                icc_sub_array[i] = icc_sub_df.loc['ICC3', 'ICC'] 

            except Warning as err:
                #print('-Error calculating. Setting as missing.')
                icc_sub_array[i] = np.NaN #If warning from ICC, turn value to missing.

            try:
                #Task-level ICC (Judges are tasks, targets are IDs)
                #Output is a dataframe
                icc_task_df = pg.intraclass_corr(data=restrict_df,
                targets='id', raters='mod', ratings='conn_value')\
                    .set_index('Type')
                icc_task_array[i] = icc_task_df.loc['ICC3', 'ICC']

            except Warning as err:
                icc_task_array[i] = np.NaN

    return icc_sub_array, icc_task_array

def reshaping(icc_sub_array, icc_task_array):
    """ Simple reshaping from 200,000x1 array to a 400x400 array
    """

    #Since order is conserved when we flatten and when we do the loop, the
    # original order is conserved when reshaping.
    res_icc_sub_array = np.reshape(icc_sub_array, (400,400))
    res_icc_task_array = np.reshape(icc_task_array, (400,400))

    return res_icc_sub_array, res_icc_task_array

def export(res_icc_sub_array, res_icc_task_array, output, name):
    """ Simple function exporting the computed arrays.
    """

    np.savetxt(f"{output}/res_icc_sub_{name}.csv", res_icc_sub_array, fmt='%1.3f', delimiter=',')
    np.savetxt(f"{output}/res_icc_task_{name}.csv", res_icc_task_array, fmt='%1.3f', delimiter=',')

if __name__ == "__main__":
	main()