#!/usr/bin/python
"""
# Fingerprinting

Script - Fingerprinting calculation
Author: Frédéric St-Onge (adapted from a script of Jordana Remz)

Date created: 2020-02-27
Date modified: 2022-10-28

--------
Purpose:
--------
    This script calculates the self- and others-identifiability, the differential
    identifiability [1] and the fingerprint identification accuracy [2]. The goal is to 
    determine whether functional connectivity matrices of an individual correlate between fMRI
    modalities. When they correlate better than with with matrices of other individuals, we
    assume that the functional connectome fingerprint is accurate.

    The script outputs a single .csv file containing the measures above. The
    others-identifiability and the differential identifiability are exported as averages.
    However, the script also outputs the similarity matrix constructed during the
    fingerprinting, which contains the correlations between the functional connectivity
    of every single participants. It also outputs a list of subjects that were used in the
    fingerprinting.

--------
Usage:
--------
    This script is a command-line script. It can be launched from the command-line by using
        python fingerprinting.py -h
    Arguments required will be indicated by -h option.

    Fundamentally, the script needs the path to two folders containing functional connectivity
    matrices and the path to where the results should be output. It also requires some string
    inputs (name of modalities, parcellation, network used, correlation type, etc.). This is
    mostly due to the needs of the original study.

    The most import of the string inputs are the parcellation and network argument as they
    will determine what brain regions (i.e., nodes) will be used, and from which parcellation.
    Currently, the script only accepts networks included in the Schaefer atlas (400 nodes) or
    the Power atlas (264 nodes). Future versions will allow more flexibility.

    A new addition to the script is the possibility to choose whether to use within- or
    between-network connections for fingeprinting. Currently, between-network is only available
    for the Schaefer atlas.

    Finally, the script allow users to specify whether only specific subjects should be
    retained for the analyses using the qc1, qc2 and reduction arguments. By default, the
    script will use all of the matrices available in the folders provided by the user. Instead,
    users can provide a list of participants (text file, 1 ID per line) of participants to 
    keep. Regex is used to match ID to file names, so names need to match exactly right.

    Particularities:
        - Due to the original study, names of matrices fed to the script need to end with 
        either ".txt" or "matrix"
        - Similarly, the script except that participant IDs will be 6 digit IDs. The script
        will not work otherwise. 

--------
Output:
--------
    The script creates a special folder architecture based on the arguments provided by the users. The goal was originally to isolate the different runs of fingerprinting from one
    another. In the output path, three folders are created; one each for the fingerprint
    metrics, the similarity matrices and the used subject list. Within each folder, subfolders
    are created, matching the different parameters used (within/between-network, network,
    parcellation, correlation type.).

--------
Future updates:
--------

    Work is currently underway to adapt this script into a Python package that would be easier
    to use and give more flexibility to users (open the choice to more atlases, offer other
    methods for fingerprinting, offer the possibility to use spreadsheets for fingerprinting,
    etc.)

--------
References:
--------

    [1]: Finn et al. (2015). Nat. Neuro.
    [2]: Amico & Goñi. (2018). Sci. Reports

"""

import os
import numpy as np
import re
import warnings
import pandas as pd
from argparse import ArgumentParser

def main():
    """
    This function launches the fingerprinting.

    """ 
    print('--------------------------------------------')
    print('---------Fingerprinting with Python---------')
    args = parsing()
    print('Parameters:')
    print(f'    Modality 1: {args.modality1} / Modality 2: {args.modality2}')
    print(f'    Path to matrices of modality 1: {args.path_m1}')
    print(f'    Path to matrices of modality 2: {args.path_m2}')
    print(f'    Path to output: {args.output}')
    print(f'    Type of fingerprinting: {args.type}')
    print(' ')
    print(f'    Parcellation used: {args.parcellation}')
    print(f'    Network used: {args.network}')
    print(f'    Correlation used: {args.corr}')
    print(f' ')
    print(f'    QC List modality 1: {args.qc_list_m1}')
    print(f'    QC List modality 2: {args.qc_list_m2}')
    print(f'    Subset: {args.reduction}')
    print(f'    Name of subset: {args.extra}')
    print('--------------------------------------------')
    print(' ')

    verify_input(args)
    print(f'Initializing an object of class fingerprint...')
    #Create an object of the Fingerprint class
    fp = fingerprint(args.path_m1, args.path_m2, args.qc_list_m1, args.qc_list_m2, args.output,
    args.type, args.modality1, args.modality2, args.parcellation, args.network, args.corr,
    args.reduction, args.extra)

    #print(f'Setting nodes of interest based on the {parcellation} parcellation, with {network} network')
    print(f'Selecting ROIs based on network and parcellation...')
    fp.rois()

    print(f'Scrubbing directories and fetching filenames...')
    fp.filename_fetch()

    print(f'Subsetting the subjects based on QC_list and reduction (if applicable)...')
    fp.subset()

    print(f'Finding out which subjects are in both modalities...')
    fp.final_subject_list()

    print(f'Targeting filenames of subjects in both modalities...')
    fp.final_file_selection()

    print(f'Calculating the similarity matrices...')
    fp.sim_matrix_calculation()

    print(f'Export the similarity matrix to a file...')
    fp.export_matrix()

    print(f'Export the subject list to a file...')
    fp.export_subject_list()

    print(f'Calculating the fingerprinting accuracy...')
    fp.fingerprint_accuracy()

    print(f'Calculating the FPC, ASC and BSD...')
    fp.fpc_asc_calculation()

    print(f'Export the FP metrics and the accuracy of the FP to a file...')
    fp.export_fpc_df()

    print('Extraction complete!')
    print('--------------------------------------------')
    print(' ')
    return None

def parsing():
    """
    Argument parser for the current script. This function:
        1. Creates a parser
        2. Add the arguments necessary to the parser object
        3. Creates an args object by parsing the arguments
        4. The args object is returned.
    """
    parser = ArgumentParser(description='Fingerprint subjects based on two different modalities. It writes out 1 similarity matrix per fingerprint and 1 subject list. It then separates the FPC from similarity coefficients and outputs them in a clean .csv file')
    parser.add_argument('--path_m1', help='Path to find connectivity matrices for modality 1')
    parser.add_argument('--path_m2', help='Path to find connectivity matrices for modality 2')
    parser.add_argument('--output', help='Path where the data should be output')
    parser.add_argument('-t', '--type', choices=['within', 'between'], default='within', help='(Optional) Type of fingerprinting (do we use within- or between-network to do the fingerprint?). By default uses within.')
    parser.add_argument('-m1', '--modality1', help='Name of first modality to be used in the fingerprinting')
    parser.add_argument('-m2', '--modality2', help='Name of second modality to be used in the fingerprinting')
    parser.add_argument('-q1', '--qc_list_m1', default='', help='(Optional) The list of subjects that passed QC in modality 1. This is generated elsewhere.')
    parser.add_argument('-q2', '--qc_list_m2', default='', help='(Optional) The list of subject that passed QC in modality 2. This is generated elsewhere.')
    parser.add_argument('-p', '--parcellation', choices=['Power', 'Schaefer'], help='Atlas used for parcellation. Currently, only Power and Schaefer atlas are supported.')
    parser.add_argument('-n', '--network', default='Whole Brain', help='Networks of interest from the chosen parcellation. By default, Whole Brain is selected.')
    parser.add_argument('-c', '--corr', help='What type of correlation was used to generate the matrices?')
    parser.add_argument('-r', '--reduction', default='', help='(Optional) Argument to import a list by which the list of subject to fingerprint is subset. The argument should be a list of IDs to keep. Note that this is done PRIOR to QC exclusions.')
    parser.add_argument('-e', '--extra', default='', help='(Optional) Extra suffix to add to the output.')

    args = parser.parse_args()

    return args
    
def verify_input(args):
    
    """
    This function verifies that the input is not problematic for the fingerprinting. If an input is wrong, a SystemExit() is raised and the program stops. Currently, the script checks for:

    - Validity of the paths given (for the inputs and output; that the paths exists and all three are different)
    - That the modality 1 is different from the modality 2
    - That the networks given as input exist in the dictionary of the script.
    """
    #List of networks from the dictionary of the fingerprint class.
    network_power = ["uncertain", "sensory_somatomotor", "cingulo_opercular", "auditory", "default_mode", "ventral_attention","visual", "frontoparietal", "salience", "cerebellar", "subcortical", "dorsal_attention", "Whole_brain"]

    network_schaefer = ["visual", "somatomotor", "dorsal_attention","salience_ventral_attention", "limbic", "frontoparietal", "default_mode", "Whole_brain", "random1", "random2"]

    # Check directories for modalities
    if args.path_m1:
        if (os.path.isdir(args.path_m1) == True):
            pass
        else:
            raise SystemExit(f'The path given for the first modality ({args.path_m1}) does not exist. Please verify the path.')
    else:
        raise SystemExit('The "--path_m1" argument is mandatory.')

    if args.path_m2:
        if (os.path.isdir(args.path_m2) == True):
            pass
        else:
            raise SystemExit(f'The path given for the second modality ({args.path_m2}) does not exist. Please verify the path.')
    else:
        raise SystemExit('The "--path_m2" argument is mandatory.')

    if args.output:
        if (os.path.isdir(args.output) == True):
            pass
        else:
            raise SystemExit(f'The path given for the output ({args.output}) does not exist. Please verify the path.')
    else:
        raise SystemExit('The "--output" argument is mandatory.')

    #Check that the directories for modality 1 and 2 are different (not the same matrices)
    if args.path_m1 == args.path_m2:
        raise SystemExit(f'The path given for both modalities is the same. Please provide different paths for each of the two modalities.')

    if ((args.path_m1 == args.output) or (args.path_m2 == args.output)):
        raise SystemExit('Please provide a different directory for the output then the input directories.')

    #Check that the modalities given are not the same for the modalities
    if args.modality1 == args.modality2:
        raise SystemExit(f'Modalities given are the same ({args.modality1} and {args.modality2}). Please give different modalities name. If the same modality is used but a different timepoint or run, please simply indicate so in the name.')
    
    #Check that the network used are usable
    if args.parcellation == "Schaefer":
        if args.network in network_schaefer:
            pass
        else:
            raise SystemExit(f'The network given for the Schaefer atlas is not recognized. Please use one among the list: {network_schaefer}')
    
    if args.parcellation == "Power":
        if args.network in network_power:
            pass
        else:
            raise SystemExit(f'The network given for the Power atlas is not recognized. Please use one among the list (except cerebellar): {network_power}')
    
    #Check that the reduction option is an existing file AND is a list of IDs
    if args.reduction:
        if (os.path.isfile(args.reduction) == True):
            pass    
        else:
            raise SystemExit(f"The file given for the list to subset the data doesn't exist or was not specified properly")

    return None

class fingerprint:
    """
    The fingerprint class is built with the arguments given by the parser.
    """
    def __init__(self, path_m1, path_m2, qc_list_m1, qc_list_m2, output, type, modality1, modality2, parcellation, network, corr, reduction, extra):
        self.path_m1 = path_m1 #Path for the first modality
        self.path_m2 = path_m2 #Path for the second modality
        self.output = output #Path for the output of the data
        self.type = type #What do we use to do the fingerprinting? Within- or between-network connections?

        self.qcm1 = qc_list_m1 #Path for the QC'ed subjects included in modality 1
        self.qcm2 = qc_list_m2 #Path for the QC'ed subjects included in modality 2

        self.m1 = modality1 #Name of the first modality
        self.m2 = modality2 #Name of the second modality
        self.parce = parcellation #Name of the parcellation used
        self.net = network #Name of the network of interest
        self.corr = corr #Correlation used to generate the similarity matrices
        self.red = reduction #List of IDs to subset by. See the helper script
        self.extra = extra #Extra suffix to add to the final name

        return None

    def rois(self):
        """
        This function defines the regions of interest to be used in this fingerprinting run.

        The function simply takes the parcellation parameter from the class to determine which dictionnary to use. Then, it uses the "network" provides to determine the nodes to use for the fingerprinting.

        WARNING: The "network" variable need to mirror exactly what is provided by the original author. Otherwise, the script will not recognize it.
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
            "cerebellar": list(range(242, 246)),
            "dorsal_attention": [250, 251] + list(range(255, 264)),
            "Whole_brain": list(range(0, 264)) #If modified Power atlas, you will need to change this back to 273
        }

        random_net_array1 = None
        random_net_array2 = None
        if self.net == "random1":
            np.random.seed(667)
            random_net_array1 = list(np.sort(np.random.randint(0,400,24)))
        if self.net == "random2":
            np.random.seed(667)
            random_net_array2 = list(np.sort(np.random.randint(0,400,91)))
            print(random_net_array2)
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
            "random1": random_net_array1,
            "random2": random_net_array2
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
        
        #If we use the "within" network fingerprinting, we only care that the nodes are exactly the same to select the matrix section
        elif self.parce == "Schaefer" and self.type == "within":
            self.nodes = dict_net_schaefer_within[f'{self.net}']
        
        #If we use the "between" network fingerprinting, we want to bind the matrix to the network of interest ("within") nodes
        # and we want to select all the nodes from that slice of network that are NOT within network.
        #So for the visual network, we want edges that are between 1 and 30 on the x-axis, but we want all edges that
        # are not in the 1-30 position in the y-axis.
        elif self.parce == "Schaefer" and self.type == "between":
            self.nodes = dict_net_schaefer_within[f'{self.net}']
            self.between = dict_net_schaefer_between[f'{self.net}']
        
        else:
            raise SystemExit(f'The combination parcellation "{self.parce}" and type "{self.type}" is not registered.')

        return None

    def filename_fetch(self):
        """
        This method searches the directories where the input is supposed to be. If the files 
        end with "matrix.txt" or with simply "matrix", we append that filename to the list of the correct modality. The filename is a string.
        """
        self.filename_list_m1 = []
        self.filename_list_m2 = []

        try:
            for filename in sorted(os.listdir(self.path_m1)): #We sort the content of the directory and list the folders
                if (filename.endswith(".txt") or filename.endswith("matrix")): #If the name matches...
                    self.filename_list_m1.append(filename) #... we append it to the filename. 
            for filename in sorted(os.listdir(self.path_m2)):
                if (filename.endswith(".txt") or filename.endswith("matrix")):
                    self.filename_list_m2.append(filename)
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

        TO DO: Could probably reduce it to a loop?
        """

        #Set the empty lists
        self.qc_list_m1 = [] #To store the IDs of modality 1 (QC'ed)
        self.qc_list_m2 = [] #To store the IDs of modality 2 (QC'ed)
        self.reduct = []     #To store the IDs of the subset to do

        #######################################
        #Import the QC_list for modality 1 (with a try/except statement to catch errors)
        if self.qcm1: #If the QC_list was selected...
            try:
                with open(self.qcm1) as x: #Open file
                    for line in x:
                        new_line = re.search(r'[0-9]{6}', line).group()[0:6] 
                        if new_line.strip().isdigit() and len(new_line.strip()) == 6: #If stripped line are only numbers and their length is 6, then...
                            self.qc_list_m1.append(new_line.strip()) #Append the number to the list
                        else:
                            raise SystemExit('Error: encountered unexpected input in QC list 1. Cannot strip to one six-digit ID (CC ID) per line.')
                print('     Generated QC_list 1')

            except OSError: #If we can't find the file, we exist straight away.
                print(f"Error: Could not access QC list file at : {self.qcm1}", flush=True)
                raise SystemExit

        else: #If no QC_list is provided, print the message
            print(f'    No argument given to -q1 argument. Assume no QC_list is necessary for modality 1')
        
        #######################################
        #Import the QC_list for modality 2 (with a try/except statement to catch errors)
        if self.qcm2:
            try:
                with open(self.qcm2) as y: #Open file
                    for line in y:
                        new_line = re.search(r'[0-9]{6}', line).group()[0:6]
                        if new_line.strip().isdigit() and len(new_line.strip()) == 6: #If stripped line are only numbers and their length is 6, then...
                            self.qc_list_m2.append(new_line.strip()) #Append the number to the list
                        else:
                            raise SystemExit('Error: encountered unexpected input in QC list 2. Cannot strip to one six-digit ID (CC ID) per line.')
                print('     Generated QC_list 2') #Will print if successful.
        
            except OSError: #If we can't find the file, we exist straight away.
                print(f"Error: Could not access QC list file at : {self.qcm2}", flush=True)
                raise SystemExit
        else:
            print(f'    No argument given to -q2 argument. Assume no QC_list is necessary for modality 2')
        
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

        (FUNCTION COULD PROBABLY BE SIMPLIFIED/MORE FLEXIBLE.)
        """
        #We extract the digits from the filename to create the subject list of m1 and m2
        self.sub_list_m1 = [re.search(r's[0-9]{6}', string).group()[1:7] for string in self.filename_list_m1]
        self.sub_list_m2 = [re.search(r's[0-9]{6}', string).group()[1:7] for string in self.filename_list_m2]

        #We set the list of IDs as the list of list of IDs from the filenames, the QC lists
        #and any subset we want to do. This is checked first.
        if self.qcm1 and self.qcm2 and self.red:
            self.all_IDs = [self.sub_list_m1, self.sub_list_m2, self.qc_list_m1, self.qc_list_m2, self.reduct]

        elif self.qcm1 and self.qcm2:
            self.all_IDs = [self.sub_list_m1, self.sub_list_m2, self.qc_list_m1, self.qc_list_m2]
        
        else:
            self.all_IDs = [self.sub_list_m1, self.sub_list_m2]

        #Find the intersection of all subset lists (equivalent of "inner" in a dataframe merge)
        self.subject_intersection_set = set(self.all_IDs[0]).intersection(*self.all_IDs)

        #Transform the dictionary of the previous variable to a sorted list
        self.subject_list = list(sorted(self.subject_intersection_set))

        return None

    def final_file_selection(self):
        """
        We create a list that will house the final filenames to fetch and use for the fingerprinting. We do this by iterating over the subject_list we created in the previous method AND iterating over the filenames available for both modalities, one at a time.

        I.e. Our subject list include subjects that are in BOTH modalities, while our list of files for each modality includes subjects ONLY in this modalities. Fingerprinting needs 2 modalities to be calculated, so we restrict the filenames to take to subjects that also have the other modality.
        """

        self.filename_list_m1_final = [filename for subject in self.subject_list for filename in self.filename_list_m1 if subject in filename]
        self.filename_list_m2_final = [filename for subject in self.subject_list for filename in self.filename_list_m2 if subject in filename]

        return None
    
    def sim_matrix_calculation(self):
        """
        This method calculates the similarity matrix for the given modalities.

        It first creates an empty matrix that has the shape of the subject_list.

        Then, iteratively, for each node of the matrix, we load the matrix of each subject, correlate it to the matrix of all other subjects, and store each of the results in the node of the similarity matrix. 

        Oct 2022: Adding the random networks caused some issues in computations. More checks were added.
        """

        self.ns = len(self.subject_list)
        self.sim_matrix = np.zeros(shape=(self.ns, self.ns))

        for i in range(self.ns):
            print(f'i={i}', flush=True)
            try:
                matrix_file_m1 = np.loadtxt(f'{self.path_m1}/{self.filename_list_m1_final[i]}', dtype=np.double)
            except ValueError:
                matrix_file_m1 = np.loadtxt(f'{self.path_m1}/{self.filename_list_m1_final[i]}', delimiter=',', dtype=np.double)

            #Set an empty variable called r1. We will use it to store the data we need for the correlation.
            #r1 = None #The empty variable is made because it will store different values depending on the type of fingerprinting done
            if self.type == 'within':
                submatrix1 = matrix_file_m1[self.nodes][:, self.nodes]
                r1 = submatrix1[np.triu_indices(len(submatrix1), k=1)] #Exclude the diagonal and retain upper triangular (within-network connections are a "similarity" graph)
                if self.net == "random2":
                    print(r1)
            elif self.type == 'between':
                submatrix1 = matrix_file_m1[self.nodes][:, self.between] #Exclude any within-network connections (retain only the rectangle of "between" connections)
                r1 = submatrix1.flatten() #The within-network, upper triangle, gives a "flattened" array. We also need to flatten our "between" network array to execute the correlation properly 
                # (otherwise the correlation between individual is exactly the same across all subjects)

            z1 = np.arctanh(r1) #Fisher normalize the remaining edges

            #If there are cells where the Fisher norm returns "Inf", replace the "Inf" by 0
            if np.count_nonzero(np.isinf(z1)) > 0:
                fin1 = np.where(np.isinf(z1), 0, z1)
            else:
                fin1 = z1.copy()

            for j in range(i, self.ns):
                try:
                    matrix_file_m2 = np.loadtxt(f'{self.path_m2}/{self.filename_list_m2_final[j]}', dtype=np.double)
                except ValueError:
                    matrix_file_m2 = np.loadtxt(f'{self.path_m2}/{self.filename_list_m2_final[j]}', delimiter=',', dtype=np.double)

                #See comments in the first loop for more info.
                #r2 = None
                if self.type == 'within':
                    submatrix2 = matrix_file_m2[self.nodes][:, self.nodes]
                    r2 = submatrix2[np.triu_indices(len(submatrix2), k=1)]
                if self.type == 'between':
                    submatrix2 = matrix_file_m2[self.nodes][:, self.between]
                    r2 = submatrix2.flatten()

                z2 = np.arctanh(r2)

                if np.count_nonzero(np.isinf(z2)) > 0:
                    fin2 = np.where(np.isinf(z2), 0, z2)
                else:
                    fin2 = z2.copy()

                #For Fisher normalized arrays, we compute the Pearson correlation between the two modalities. This results in a similarity matrix where
                #   the diagonal is the "within-individual" correlation while off-diagonal elements are the "between-individual" correlation
                self.sim_matrix[i, j] = np.corrcoef(fin1, fin2)[0, 1] #In the future, will need to change to scipy for this. Otherwise, can't use other correlation methods for FP.
        
        #The sim matrix is computed for the upper triangle only for the within-network
        # We add a final transpose to make a full matrix, as the similarity matrix is
        # completly symmetric.
        self.sim_matrix = self.sim_matrix + np.triu(self.sim_matrix, k=1).T

        return None
    
    def export_matrix(self):
        """
        This function exports the similarity matrix. In the output, it creates a special folder (if it does not exist) where the output is stored. 
        """
        self.matrix_filename = f'similarity_matrix_{self.m1}_{self.m2}_{self.type}_{self.parce}_{self.net}_{self.corr}'

        self.path_matrix_final = f'{self.output}/similarity_matrices/{self.type}/{self.parce}/{self.corr}'
        if not os.path.exists(self.path_matrix_final):
            os.makedirs(self.path_matrix_final)

        with open(f'{self.path_matrix_final}/{self.matrix_filename}', 'w') as f:
            np.savetxt(f, self.sim_matrix, fmt='%1.3f')
            f.flush
    
        return None

    def export_subject_list(self):
        """
        This method exports the final subject list used in the fingerprinting. It creates a special folder in the output path to generate this. 
        """
        self.subject_filename = f'subject_list_{self.m1}_{self.m2}_{self.type}_{self.parce}_{self.net}_{self.corr}'

        self.path_subject_final = f'{self.output}/subject_list/{self.type}/{self.parce}/{self.corr}'
        if not os.path.exists(self.path_subject_final):
            os.makedirs(self.path_subject_final)
        
        with open(f'{self.path_subject_final}/{self.subject_filename}', 'w') as f:
            f.write('\n'.join(self.subject_list))
            f.flush()
        
        return None

    def fingerprint_accuracy(self):
        """
        This function tests which fingerprint are successful. It does this by checking if the diagonal of the matrix is the maximum (so if the correlation is strongest with self. If so, the subject is fingerprinted.). We compare the list of fingerprinted subjects to the original list. This generates a list of subjects that were not fingerprinted, which we use to build a new dataframe. 
        """
        fingerprinted_list = []

        #For every subject in our dataset...
        for i in range(self.ns):
            #If the maximum value is achieved at the diagonal...
            if np.argmax(self.sim_matrix[i, :]) == i and np.argmax(self.sim_matrix[i, :][::-1]) == (self.ns - 1 - i):
                #Append the subject ID to a "fingerprinted" subject_list.
                fingerprinted_list.append(self.subject_list[i])
        
        #Create a dataframe with 2 columns: the index AND the FP status.
        non_fp_dataframe = pd.DataFrame(columns=[f'participant_id', f'non_fp_{self.m1}_{self.m2}_{self.type}_{self.parce}_{self.net}_{self.corr}'])

        #Use a condition to check if, for each row, the ID is in the fingerprinted list. If the subject is fingerprinted, we assign a 1. Otherwise, assign a 0.
        non_fp_dataframe['participant_id'] = self.subject_list
        non_fp_dataframe[f'non_fp_{self.m1}_{self.m2}_{self.type}_{self.parce}_{self.net}_{self.corr}'] = np.where(non_fp_dataframe['participant_id'].astype("str").isin(fingerprinted_list), 1, 0)

        #Set a final index to the dataframe, which will be "participant_id"
        self.non_fp_dataframe = non_fp_dataframe.set_index('participant_id')

        return None

    def fpc_asc_calculation(self):
        """
        This function is to create the dataframes containing three metrics: the Fingerprinting
        Coefficient (FPC), the Average Similarity Coefficient (ASC) and the Between-Subject Distance (BSD).

        Note that this naming is old and was never changed in the code. However, Fingerprinting
        Coefficient (FPC) corresponds to self-identifiability, average similarity coefficient
        (ASC) corresponds to others-identifiability and between-subject distance (BSD)
        corresponds to differential identifiability. New scripts 

        The FPC is equivalent to the correlation "within-subject". It is extracted by
        extracting the diagonal of the similarity matrix.
        The ASC is equivalent to the average "between-subject" correlation. It is extracted by
        averaging all the cells in a single row, minus the diagonal cell (within-subject).
        The BSD is the Eucledian Distance between the FPC and the ASC.

        The function:
        1. Calculates all three metrics and saves them as instances of our object
        2. Creates a dictionnary which will save the three columns and their values (list to column)
        3. Creates a final dataframe where the index is the subject_list, and the data is the
        dictionnary (so all three coefficients with correct column names)
        """ 
        #We calculate the metrics to be used
        self.fpc_coef = np.diag(self.sim_matrix) #Within-individual correlation
        self.asc_coef = (self.sim_matrix.sum(1)-np.diag(self.sim_matrix))/(self.sim_matrix.shape[1]-1) #Average between individual correlation
        self.identi = self.fpc_coef - self.asc_coef #Identifiability (based on Amico 2018), which is simply fpc - asc.

        #We create a dictionary we will use to populate our dataframe
        self.coef_dict = {f'fpc_{self.m1}_{self.m2}_{self.type}_{self.parce}_{self.net}_{self.corr}':self.fpc_coef,
        f'asc_{self.m1}_{self.m2}_{self.type}_{self.parce}_{self.net}_{self.corr}':self.asc_coef,
        f'bsd_{self.m1}_{self.m2}_{self.type}_{self.parce}_{self.net}_{self.corr}':self.identi} 
        #Abreviations in order: FPC = Fingerprinting coefficient, ASC =  Average Similarity Coefficient
        #BSD = Between subject_distance

        #We create the dataframe, using the dictionnary for data/columns and the list of subjects 
        #as the index
        self.sim_matrix_dataframe = pd.DataFrame(self.coef_dict, index=self.subject_list)

        return None 
    
    def export_fpc_df(self):
        """
        This final function merges the FPC dataframe and the fingerprinting accuracy dataframe to result in a single dataframe.
        """

        #Merge the calculated FPC/ASC/Iden values and merge them to the fingerprint status
        self.final_fpc_df = self.sim_matrix_dataframe.merge(self.non_fp_dataframe, left_index=True, right_index=True, how='left')

        #Sets a name for the index during export
        self.final_fpc_df.index.name = f'participant_id'

        #If we need to give an extra name to the file, it is added here. Otherwise, a name without extra is outputed.
        if self.extra:
            self.final_fpc_filename = f'fpc_{self.m1}_{self.m2}_{self.type}_{self.parce}_{self.net}_{self.corr}_{self.extra}'
        else:
            self.final_fpc_filename = f'fpc_{self.m1}_{self.m2}_{self.type}_{self.parce}_{self.net}_{self.corr}'

        #Check if the output path exists. If not, create it.
        path_fpc_final = f'{self.output}/fpc/{self.type}/{self.parce}/{self.corr}'
        if not os.path.exists(path_fpc_final):
            os.makedirs(path_fpc_final)
        
        self.final_fpc_df.to_csv(f'{path_fpc_final}/{self.final_fpc_filename}.csv')

        return None
		
	
if __name__ == "__main__":
	main()
