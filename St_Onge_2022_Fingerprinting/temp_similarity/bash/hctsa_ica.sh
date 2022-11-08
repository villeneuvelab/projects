#!/bin/bash

#Set the pairs of modality we want to run with the R script.
declare -a StringArray=("--modality rest" "--modality task" "--modality movie")

echo "--------------------------"
echo "R Script - Var with ICAs"
echo "--------------------------"

#We run the R script for each pair of modality using a for loop
for STRING in "${StringArray[@]}"; 
do 
    echo "--------------------------"
    echo "Processing pair ${STRING}"
    echo " "
    
    Rscript var_ica.R --path_file /Users/stong3/Desktop/Fingerprinting_Project/CAMCAN/05_Stats/01_Datasets/fpc_schaefer_pcor_data_clean.csv --path_ica /Users/stong3/Desktop/Fingerprinting_Project/CAMCAN/05_Stats/01_Datasets/merged_all_GM_CamCAN_ica_common_template_no_emci_Demographics.csv --output /Users/stong3/Desktop/Fingerprinting_Project/CAMCAN/05_Stats/03_New_stats/hctsa_ica --atlas Schaefer --corr hctsa $STRING

done