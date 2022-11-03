#!/bin/bash

#Set the pairs of modality we want to run with the R script.
#"--modality1 rest --modality2 movie" "--modality1 task --modality2 movie" "--modality1 task --modality2 rest" "--modality1 movie --modality2 rest" "--modality1 movie --modality2 task"
declare -a Modalities=("--modality1 rest --modality2 task")
declare -a Type=("--type within" "--type between")

echo "--------------------------"
echo "R Script - FP Identifiability - Schaefer"
echo "--------------------------"

#We run the R script for each pair of modality using a for loop
for STRING in "${Modalities[@]}"; do
    echo "--------------------------"
    echo "Processing pair ${STRING}"
    echo " "
    for TYPE in "${Type[@]}"; do

        echo "For type ${TYPE}"
        Rscript fpc_iden.R --path_file /ABSOLUTE/PATH/TO/FINGERPRINT/DATA/fp_data_Schaefer.csv  --output /ABSOLUTE/PATH/TO/OUTPUT --atlas Schaefer --corr pcor $STRING $TYPE
    done

done

echo "--------------------------"
echo "R Script - FP Identifiability - Power"
echo "--------------------------"

#We run the R script for each pair of modality using a for loop
for STRING in "${Modalities[@]}"; do
    echo "--------------------------"
    echo "Processing pair ${STRING}"
    echo " "
    for TYPE in "${Type[@]}"; do

        echo "For type ${TYPE}"
    
        Rscript fpc_iden.R --path_file /ABSOLUTE/PATH/TO/FINGERPRINT/DATA/fp_data_Schaefer.csv  --output /ABSOLUTE/PATH/TO/OUTPUT --atlas Power --corr pcor $STRING $TYPE
    done

done