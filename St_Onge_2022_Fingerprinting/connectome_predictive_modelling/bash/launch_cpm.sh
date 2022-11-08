#!/bin/bash

#Simple bash file to launch the CPM.
# The goal is to be able to launch the script and let it run, without needing anything from my end
# Due to the nature of the time constraints, I am doing this a bit quick and dirty

#The script below assumes that the subject lists are stored and labeled identically. 
#The script also has no "checkpoint" per say. Once something is output, it goes on to the next CPM.
#The python script doesn't check if the data exists before running.

#Add more to the lists to run more parameters, modalities with one call
declare -a modality=("rest")
declare -a win_size=("100")
declare -a step_size=("40")

#Other variables to help with the script writting
PYT_SCRIPT="/PATH/TO/PYTHON/SCRIPT/cpm_python_v2.py"
OUTPUT="/PATH/TO/OUTPUT"

#Loop over the lists to launch the different CPMs
for mod in ${modality[@]}
do
    for win in ${win_size[@]}
    do
        for ss in ${step_size[@]}
        do
        COUNTER=0
            for FILE in /PATH/TO/WINDOW/DATA/"${win}"_"${ss}"/*
            do
            let COUNTER++
            echo "##########"
            if (("${COUNTER}" < 10)); then
                echo "Working on win 0${COUNTER} for modality ${mod} with size ${win} and step ${ss}"
                python ${PYT_SCRIPT} -p1="/PATH/TO/MATRICES/adj_matrix_Schaefer_${mod}_n400_pcor" -d="${FILE}" -o="${OUTPUT}"/"${win}_${ss}" -n="fpc_${mod}_${win}_${ss}_win0${COUNTER}" -m=LM
                echo "Done win 0${COUNTER}!"
            else
                echo "Working on win ${COUNTER} for modality ${mod} with size ${win} and step ${ss}"
                python ${PYT_SCRIPT} -p1="/PATH/TO/MATRICES/adj_matrix_Schaefer_${mod}_n400_pcor" -d="${FILE}" -o="${OUTPUT}"/"${win}_${ss}" -n="fpc_${mod}_${win}_${ss}_win${COUNTER}" -m=LM
                echo "Done win ${COUNTER}!"
            fi

            echo "#########"

            done
        done
    done
done
