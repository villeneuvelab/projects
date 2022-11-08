"""
# Fingerprinting - CAMCAN Cohort
## Variance in functional connectivity

Script - Variance calculator
Author: Frédéric St-Onge

Date created: 2021-07-23
Date modified: 2022-01-25

Version: 1.1

--------
Updates:
--------
    2021-12-14: Added documentation across. Added the calculation of absolute distance from mean (in addition to the median stuff) so we can look at both.
    2022-01-25: Added the possibility to look at variability "between" network as well (not just within)

--------
Purpose:
--------
    This script computes the variance and the coefficient of variance of the functional connectivity matrices of each subject.
    It also computes the distance between the individual-level variance and the group-level mean and median of individual-level variance.
    It does so for each network of interest.

    This script will work for any symmetric matrix matching what we would expect for functional connectivity. As such, it also works
    for the temporal similarity analysis.

"""

import os
import numpy as np
import re
import warnings
import pandas as pd
from argparse import ArgumentParser

def main():
    """
    This function launches the variance_calculator.


    """ 
    args = parsing()
    print('--------------------------------------------')
    print('----Connectivity variability with Python----')
    print('Parameters:')
    print(f'    Path to matrices, modality {args.modality}: {args.path_mat}')
    print(f'    Path to output: {args.output}')
    print(f'    Index name: {args.index}')
    print(f'    Type of network used: {args.type}')
    print(' ')
    print(f'    Modality of interest: {args.modality}')
    print(f'    Parcellation used: {args.parcellation}')
    print(f'    Network used: {args.network}')
    print(f'    Correlation used to generate the original matrix: {args.corr}')
    print(f'    QC List used: {args.qcm}')
    print(f'    Additional participant selection used: {args.red}')
    print('--------------------------------------------')
    print(' ')

    #verify_input(args)
    print(f'Initializing an object of class var_calculator...')
    var = var_calculator(args.path_mat, args.output, args.type, args.modality, args.parcellation, args.network, args.corr, args.index, args.qcm, args.red)

    print(f'Selecting ROIs')
    var.rois()

    print(f'Fetching file names')
    var.filename_fetch()

    print(f'Calculating variance and exporting to file')
    var.var_calculation()

    print(f'Computing the distance between variance coefficient and the group-level median')
    var.absolute_distance()

    print(f'Exporting variance to file')
    var.export_df()

    print("Done!")

    return None

def parsing():
    """
    Argument parser for the current script.

    """
    parser = ArgumentParser(description='')
    parser.add_argument('--path_mat', help='Path to matrices to analyse')
    parser.add_argument('--output', help='Path where the output should go')
    parser.add_argument('-t', '--type', default = 'within', choices=['within', 'between'], help='Whether we should use within- or between-network edges.')
    parser.add_argument('-m', '--modality', help='Name of the modality to use. Purely used to add to the final name.')
    parser.add_argument('-p', '--parcellation', choices=['Power', 'Schaefer'], help='Which parcellation to use.')
    parser.add_argument('-n', '--network', default='Whole Brain', help='Which network to use')
    parser.add_argument('-c', '--corr', choices=['corr', 'pcor', 'hctsa'], default='corr', help='String to add, representing the correlation used. It also adds a string for "hctsa" if the script is used that way.')
    parser.add_argument('-i', '--index', default='id_camcan', help='This string should be the name of the column you want to set as index name for the resulting dataframe. Currently, by default, it will be "id_camcan".')
    parser.add_argument('-q', '--qcm', default='', help="List of QC'ed participants (participants that should be retained).")
    parser.add_argument('-r', '--red', default='', help='Further selection of subject, if needed.')

    args = parser.parse_args()

    return args
    

class var_calculator:

    def __init__(self, path_mat, output, type, modality, parcellation, network, corr, index, qcm, red):
        self.path_mat = path_mat
        self.output = output
        self.type = type

        self.modality = modality
        self.parce = parcellation
        self.net = network
        self.corr = corr
        self.index = index
        self.qcm = qcm
        self.red = red

        return None

    def rois(self):
        """
        This function defines the regions of interest to be used in this connectivity variability run.

        Currently, the script supports the Power (REF) and the Schaefer atlas (REF). Note that the Power atlas in our cohort have also been processed with an additional 8 nodes (limbic), based on a paper by Vachon-Presseau et al. (REF).

        The function simply takes the parcellation parameter from the class to determine which dictionnary to use. Then, it uses the "network" to determine the nodes to use for the calculation.

        WARNING: The "network" variable need to mirror exactly what is written below. Otherwise, the script will not recognize it.

        FEATURE TO ADD: Create a way where we can add different networks together. For now, any "combos" need to be added in the dictionnary
        """
        #Dictionnary for the Power Atlas
        dict_net_power = {
            #REDO INDEXING. FIRST TERM SHOULD BE -1
            "uncertain": list(range(0, 12)) + [83, 84] + list(range(131, 136)) + list(range(139, 142)) + list(range(181, 185)) + [220] + list(range(246, 250)) + [252, 253],
            "sensory_somatomotor": list(range(12, 46)) + [254],
            "cingulo_opercular": list(range(46, 60)),
            "auditory": list(range(60, 73)),
            "default_mode": list(range(73, 83)) + list(range(85, 131)) + [136] + [138],
            "ventral_attention": [137] + list(range(234, 242)),
            "visual": list(range(142, 173)),
            "frontoparietal": list(range(173, 181)) + list(range(185, 202)),
            "salience": list(range(202, 220)),
            "subcortical": list(range(221, 234)),
            #"cerebellar": list(range(242, 246)),
            "dorsal_attention": [250, 251] + list(range(255, 264)),
            "limbic": list(range(264, 273)),
            "default_mode_limbic": list(range(73, 83)) + list(range(85, 131)) + [136] + [138] + list(range(264, 273)),
            "Whole_brain": list(range(0, 273))
        }

        #Dictionnary for the Schaefer Atlas (Within- network)
        dict_net_schaefer_within = {
            #Based on the dataset.labels from the Schaefer Atlas in Nilearn.
            "visual": list(range(0, 31)) + list(range(200, 230)),
            "somatomotor": list(range(31, 68)) + list(range(230, 270)),
            "dorsal_attention": list(range(68, 91)) + list(range(270, 293)),
            "salience_ventral_attention": list(range(91, 113)) + list(range(293, 318)),
            "limbic": list(range(113, 126)) + list(range(318, 331)),
            "frontoparietal": list(range(126, 148)) + list(range(331, 361)),
            "default_mode": list(range(148, 200)) + list(range(361, 400)),
            "Whole_brain": list(range(0, 400)),
            "default_mode_limbic": list(range(148, 200)) + list(range(361, 400)) + list(range(113, 126)) + list(range(318, 331))
        }

        #Dictionary for the Schaefer Atlas (Between- network) 
        dict_net_schaefer_between = {
            #Overall between-network connectivity (we don't do the average between each network)
            "visual": list(range(31, 200)) + list(range(231, 400)),
            "somatomotor": list(range(0, 31)) + list(range(68, 230)) + list(range(270, 400)),
            "dorsal_attention": list(range(0, 68)) + list(range(91, 270)) + list(range(293, 400)),
            "salience_ventral_attention": list(range(0, 91)) + list(range(113, 293)) + list(range(318, 400)),
            "limbic": list(range(0, 113)) + list(range(126, 318)) + list(range(331, 400)),
            "frontoparietal": list(range(0, 126)) + list(range(148, 331)) + list(range(361, 400)),
            "default_mode": list(range(0, 148)) + list(range(200, 361))
            }
        if self.parce == "Power":
            self.nodes = dict_net_power[f'{self.net}']
        
        #If we use the "within" network FPC, we only care that the nodes are exactly the same to select the matrix section
        elif self.parce == "Schaefer" and self.type == "within":
            self.nodes = dict_net_schaefer_within[f'{self.net}']
        
        #If we use the "between" network FPC, we want to bind the matrix to the network of interest ("within") nodes
        # and we want to select all the nodes from that slice of network that are NOT within network.
        #So for the visual network, we want edges that are between 1 and 30 on the x-axis, but we want all edges that
        # are not in the 1-30 position in the y-axis.
        elif self.parce == "Schaefer" and self.type == "between":
            self.nodes = dict_net_schaefer_within[f'{self.net}']
            self.between = dict_net_schaefer_between[f'{self.net}']
        
        return None

    def filename_fetch(self):
        """
        This method searches the directories where the input is supposed to be. If the files end
        with "matrix.txt" (for the Power Atlas) or with simply "matrix" (for the Schaefer Atlas),
        we append that filename to the list of the correct modality. 
        """
        self.filename_list = []

        try:
            for filename in sorted(os.listdir(self.path_mat)): #We sort the content of the directory and list the folders
                if (filename.endswith("matrix.txt") or filename.endswith("matrix")): #If the name matches...
                    self.filename_list.append(filename) #... we append it to the filename. 
        except OSError:
            print("Error: Could not access data directory.", flush=True)
            raise SystemExit
        
        return None

    def subset(self):
        """
        This method imports the files that we need to subset from the fingerprinting list.
        This is because, in our study, the connectivity matrix of participants is computed
        whether they passed the QC or not. 

        If the folder containing your connectivity matrices only holds participants who
        have been QC'ed, the -q1 and -q2 options should be left empty. This step is then
        skipped.

        If you wish to subset the fingerprinting in a specific way (e.g. male/female), and
        you have entered a list in the -r option, the fetching will also occur here.
        """

        #Set the empty lists
        self.qc_list = [] #To store the IDs of modality 1 (QC'ed)
        self.reduct = []     #To store the IDs of the subset to do

        #######################################
        #Import the QC_list for modality 1 (with a try/except statement to catch errors)
        if self.qcm: #If the QC_list was selected...
            try:
                with open(self.qcm) as x: #Open file
                    for line in x:
                        new_line = re.search(r'[0-9]{6}', line).group()[0:6] 
                        if new_line.strip().isdigit() and len(new_line.strip()) == 6: #If stripped line are only numbers and their length is 6, then...
                            self.qc_list.append(new_line.strip()) #Append the number to the list
                        else:
                            raise SystemExit('Error: encountered unexpected input in QC list. Cannot strip to one six-digit ID (CC ID) per line.')
                print('     Generated QC_list')

            except OSError: #If we can't find the file, we exist straight away.
                print(f"Error: Could not access QC list file at : {self.qcm1}", flush=True)
                raise SystemExit

        else: #If no QC_list is provided, print the message
            print(f'    No argument given to -q argument. Assume no QC_list is necessary')
        
        #######################################
        #Subset the data by a specific criteria
        if self.red: #If the reduction option is given
            try:
                with open(self.red) as z: #Open file
                    for line in z:
                        new_line = re.search(r'[0-9]{6}', line).group()[0:6]
                        if new_line.strip().isdigit() and len(new_line.strip()) == 6: #If stripped line are only numbers and their length is 6, then...
                            self.reduct.append(new_line.strip()) #Append the number to the list
                        else:
                            raise SystemExit('Error: encountered unexpected input in QC list 1. Expected one six-digit ID (CC ID) per line.')

                print('     Generated subset list')

            except OSError: #If we can't find the file, we exist straight away.
                print(f"Error: Could not access subset list file at : {self.red}", flush=True)
                raise SystemExit
        else:
            print(f"     No argument given to -r argument. Assume we don't want to subset the data and use all subjects that were QC'ed.")
        

        return None


    def final_subject_list(self):
        """
        We define the final subject list as the intersection of subjects in:
        - The filenames of the matrices of modality 1
        - The filenames of the matrices of modality 2
        - The subject IDs of the QC'ed subjects in modality 1
        - The subject IDs of the QC'ed subjects in modality 2
        Optional:
            - If we want to subset, we will do it here.

        """
        #We extract the digits from the filename to create the subject list of m1 and m2
        self.sub_list = [re.search(r's[0-9]{6}', string).group()[1:7] for string in self.filename_list]

        #We set the list of IDs as the list of list of IDs from the filenames, the QC lists
        #and any subset we want to do. This is checked first.
        if self.qcm and self.red:
            self.all_IDs = [self.sub_list, self.qc_list, self.reduct]

        elif self.qcm:
            self.all_IDs = [self.sub_list, self.qc_list]
        
        else:
            self.all_IDs = self.sub_list

        #Find the intersection of all subset lists (equivalent of "inner" in a dataframe merge)
        self.subject_intersection_set = set(self.all_IDs[0]).intersection(*self.all_IDs)

        #Transform the dictionary of the previous variable to a sorted list
        self.subject_list = list(sorted(self.subject_intersection_set))

        return None

    def final_file_selection(self):
        """
        We create a list that will house the final filenames to fetch and use for the
        fingerprinting. We do this by iterating over the subject_list we created in
        the previous method AND iterating over the filenames available for both
        modalities, one at a time.

        I.e. Our subject list include subjects that are in BOTH modalities, while our
        list of files for each modality includes subjects ONLY in this modalities. 
        Fingerprinting needs 2 modalities to be calculated, so we restrict the filenames
        to take to subjects that also have the other modality.
        """

        self.filename_list_final = [filename for subject in self.subject_list for filename in self.filename_list if subject in filename]

        return None

    def var_calculation(self):
        """
        This method calculates the variance coefficient for each subject and returns it within a dataframe.
        We compute different measures: a simple variance and a coefficient of variance. The simple variance
        is reported in the paper.
        """

        #We create an empty dataframe to store the results
        self.var_df = pd.DataFrame(columns=['id_camcan', f'var_{self.modality}_{self.type}_{self.parce}_{self.net}_{self.corr}', f'coef_var_ind_{self.modality}_{self.type}_{self.parce}_{self.net}_{self.corr}'])

        for i in range(len(self.filename_list)):
            print(f'i={i}', flush=True)

            #Extracting subject ID from filenames
            sub_id = re.search(r'[0-9]{6}', self.filename_list_final[i]).group()[0:6]
            
            #Loading matrix, selecting nodes and selecting lower triangular
            matrix_file = np.loadtxt(f'{self.path_mat}/{self.filename_list[i]}', dtype=np.double)

            #Slice the matrix to the appropriate location
            r1 = None
            if self.type == 'within':
                submatrix1 = matrix_file[self.nodes][:, self.nodes]
                r1 = submatrix1[np.triu_indices(len(submatrix1), k=1)]
            elif self.type == 'between':
                submatrix1 = matrix_file[self.nodes][:, self.between]
                r1 = submatrix1.flatten()

            #Calculating the variance coefficient
            var_coef = np.var(r1)

            #Calculating the coefficient of variation (at the individual level)
            coef_var = np.std(r1) / np.mean(r1)

            #Setting the variables as a dictionary and append the row to the dataframe
            new_row = {f'{self.index}': sub_id, f'var_{self.modality}_{self.type}_{self.parce}_{self.net}_{self.corr}': var_coef, f'coef_var_ind_{self.modality}_{self.type}_{self.parce}_{self.net}_{self.corr}': coef_var}

            self.var_df = self.var_df.append(new_row, ignore_index=True)

        return None
    
    def absolute_distance(self):
        """We try to come up with a "between" individual measure of the variance coefficient.

        Let's try simply using the distance (as would be implemented in a Mean Average Distance calculation).
        We compute both the mean and median.
        """

        #First step, get the average of the variable of interest (i.e., the coefficient of variation at the individual level)
        #avg_coef_var_ind = self.var_df[f'coef_var_ind_{self.modality}_{self.parce}_{self.net}_{self.corr}'].mean()
        med_coef_var_ind = self.var_df[f'var_{self.modality}_{self.type}_{self.parce}_{self.net}_{self.corr}'].median()
        mean_coef_var_ind = self.var_df[f'var_{self.modality}_{self.type}_{self.parce}_{self.net}_{self.corr}'].mean()

        #Next, we create a column where we will store the "absolute" distance of each individual from the median.
        # In Pandas, this is simply the substraction of a Series by a scalar, which results in a Series, which we force to absolute values. 
        self.var_df[f'med_dist_{self.modality}_{self.type}_{self.parce}_{self.net}_{self.corr}'] = (self.var_df[f'var_{self.modality}_{self.type}_{self.parce}_{self.net}_{self.corr}'] - med_coef_var_ind).abs()
        self.var_df[f'mean_dist_{self.modality}_{self.type}_{self.parce}_{self.net}_{self.corr}'] = (self.var_df[f'var_{self.modality}_{self.type}_{self.parce}_{self.net}_{self.corr}'] - mean_coef_var_ind).abs()
        #We add this as a column of the dataframe computed in the previous step. That way, no need to do any other modifications downstream.

        return None
    
    def export_df(self):
        """
        This function simply exports the connectivity variability dataframe in a single .csv file. The script creates a folder where to store the output.
        """

        #Creating filename
        final_var_filename = f'variability_{self.modality}_{self.type}_{self.parce}_{self.net}_{self.corr}'

        final_df = self.var_df.set_index(f"{self.index}")

        #Verifying and creating path to output
        path_var_final = f'{self.output}/{self.type}/{self.parce}/{self.corr}'
        if not os.path.exists(path_var_final):
            os.makedirs(path_var_final)
        
        #Output the csv to file
        final_df.to_csv(f'{path_var_final}/{final_var_filename}.csv')

        return None

if __name__ == "__main__":
	main()