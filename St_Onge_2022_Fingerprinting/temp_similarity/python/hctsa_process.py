"""
# Fingerprinting - CAMCAN Cohort
## HCTSA Similarity variance

Script - hctsa similarity
Author: Frédéric St-Onge

Date created: 2021-12-08
Date modified: 

Version: 1.0

--------
Purpose:
--------
    This script is made to compute temporal similarity like in the paper by Golia et al. (2020; eLife). In this paper,
    they investigate how BOLD signal variables are similar between brain regions. This method ressembles functional connectivity calculation,
    but uses timeseries features rather than syncronisation of signal.

    In our case, what we want from this script is to extract the similarity in BOLD signal features between each pair of regions, for each individual.
    The file that is outputted by hctsa is a NxK matrix, where N is the number of brain regions and K is the number of features.

"""

import os
import numpy as np
import re
import pandas as pd
from argparse import ArgumentParser

def main():
    """
    This function launches the HCTSA process scirpt. 


    """ 
    args = parsing()
    print('--------------------------------------------')
    print('----------Processing HCTSA output-----------')
    print('Parameters:')
    print(f'    Path to hctsa features: {args.path_vec}')
    print(f'    Modality: {args.modality}')
    print(f'    Path to output: {args.output}')
    print(f'    Regex provided by user: {args.regex_in}')
    print(' ')
    print(f'    Modality of interest: {args.modality}')
    print(f'    Parcellation used: {args.parcellation}')
    print('--------------------------------------------')
    print(' ')

    #verify_input(args)
    print(f'Initializing an object of class var_calculator...')
    hctsa = hctsa_out(args.path_vec, args.output, args.modality, args.parcellation, args.regex_in)

    print(f'Fetching file names')
    hctsa.filename_fetch()

    print(f'Calculating variance and exporting to file')
    hctsa.temp_prof_sim()

    print(f'Exporting variance to file')
    hctsa.export_mat()

    return None

def parsing():
    """
    Argument parser for the current script. This function:

    """
    parser = ArgumentParser(description='')
    parser.add_argument('--path_vec', help='Path to the vectors output by HCTSA.')
    parser.add_argument('--output', help='Path where matrices should be output.')
    parser.add_argument('-m', '--modality', help='Which modality to use (purely for name of the file at the end.)')
    parser.add_argument('-p', '--parcellation', choices=['Power', 'Schaefer'], help='Which parcellation to use (Right now only Power and Schaefer.')
    parser.add_argument('-r', "--regex_in", default="[0-9]{6}", help='Regex expression defined by the user to grab the participant. The script will group the strings of each character together.')

    args = parser.parse_args()

    return args

class hctsa_out:

    def __init__(self, path_vec, output, modality, parcellation, reg):
        """Creates an object of class `hctsa_out`

        Args:
            path_vec (str): Path to the hctsa feature-represented timeseries
            output (str): Path where the output should go to (ideally different than the input path)
            modality (str): String indicating which fMRI modality is used. Mostly for identification purposes.
            parcellation (str): String indicating which parcellation to use

        Returns:
            None: Assigning this to a variable in the main function will create the object.
        """
        self.path_vec = path_vec
        self.output = output

        self.modality = modality
        self.parce = parcellation
        self.regex = reg

        return None

    def filename_fetch(self):
        """
        This method searches the directories where the input is supposed to be. If the files end with "matrix.txt"
        (for the Power Atlas) or with simply "matrix" (for the Schaefer Atlas), we append that filename to the list
        of the correct modality. 
        """
        self.filename_list = []

        try:
            for filename in sorted(os.listdir(self.path_vec)): #We sort the content of the directory and list the folders
                if (filename.endswith(".csv")): #If the name matches...
                    self.filename_list.append(filename) #... we append it to the filename. 
        except OSError:
            print("Error: Could not access data directory.", flush=True)
            raise SystemExit
        
        return None

    def temp_prof_sim(self):
        """This method generates the timeseries inter-regional temporal profile similarities (based on Golia et al. (2020; eLife)).

        After creating an empty dictionary, the function iterates over all the files of the subjects and:
        1. Extracts the ID from the file name
        2. Loads the appropriate hctsa vector
        3. Computes row-wise Pearson correlation for a participant (correlation of timeseries features)
        4. Check the validity of the output
        5. Saves the output to a dictionnary

        """

        #We create an empty dictionnary to store the results
        self.dict_sim_matrices = {}

        #We loop over each participant file
        for i in range(len(self.filename_list)):

            #Extracting subject ID from filenames
            user_regex_string = re.compile(self.regex) #This takes the string given by the user and compiles it (makes it possible to use as an 'r' string)
            try:
                sub_id = re.search(user_regex_string, self.filename_list[i]).group() #Finds the participant ID based on the regex provided. Groups the output to a single string.
            except AttributeError: #If the regex doesn't match to anything, it throws an AttributeError (NoneType has no "group" function).
                print(f'Regex error: Could not find participant from file {self.filename_list[i]} matching to the regex {self.regex}') #We print the file that didn't work, and reprint the regex
                continue #We go over to the next loop if this error occurs

            #Print so that the user can see the cases processed
            print(f'Operation={i} / Computing subject {sub_id}', flush=True)
            
            #Loading vector
            try:
                matrix_file = np.loadtxt(f'{self.path_vec}/{self.filename_list[i]}', delimiter=",", dtype=np.double, encoding='utf-8-sig') #I added the encoding because it caused issues locally.
            except ValueError: #This was put in place when I tested correlation arrays with missing value
                print(f"The script does not support vectors where strings or missing values are present. Participant {sub_id} will not be extracted")
                continue #We go over to the next loop if this error occurs      

            #Calculating the correlation coefficient (similarity) matrix for the participant
            feat_sim = np.corrcoef(matrix_file)

            #Check that the resulting matrix is the right size.
            assert (matrix_file.shape[0] == feat_sim.shape[0]), f"Resulting correlation matrix of participant {sub_id} has a different number of rows (brain regions) than original array"
            assert (matrix_file.shape[0] == feat_sim.shape[1]), f"Resulting correlation matrix of participant {sub_id} has a different number of rows (brain regions) than expected (number of rows in original array)"

            #Create a new dictionary entry, where the key is the subject ID and the value is the correlation matrix
            self.dict_sim_matrices[f'{sub_id}'] = feat_sim

            #Print that the loop is done for a given subject.
            print(f'Operation={i} / Subject {sub_id} done.')

        return None
    
    def export_mat(self):
        """This function exports the similarity matrices of each participant to a folder.
        """

        #First, create a folder to store the output
        path_sim_final = f'{self.output}/{self.parce}'
        if not os.path.exists(path_sim_final):
            os.makedirs(path_sim_final)
        
        #Second, iterate over the dictionary created in the temp_prof_sim method and export each matrix
        for sub, matrix in self.dict_sim_matrices.items():

            #Create the filename
            filename_sub = f"{sub}_{self.modality}_hctsa_sim_matrix.txt"

            #Save the file of the participant to a file
            np.savetxt(f"{path_sim_final}/{filename_sub}", matrix, fmt='%1.3f') #We force the output to max of 3 digits.

        return None

if __name__ == "__main__":
	main()