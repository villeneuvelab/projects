#!/bin/bash
PATH_CONN="/PATH/TO/CONNECTIVITY/MATRICES"
PATH_QC="/PATH/TO/SUBJECT/LISTS"

parallel echo ::: "-m1=rest -m2=task --path_m1=${PATH_CONN}/adj_matrix_Schaefer_rest_n400_pcor --path_m2=${PATH_CONN}/adj_matrix_Schaefer_task_n400_pcor -q1=${PATH_QC}/final_list_subset_qc_rest -q2=${PATH_QC}/final_list_subset_qc_task" 
"-m1=task -m2=rest --path_m1=${PATH_CONN}/adj_matrix_Schaefer_task_n400_pcor --path_m2=${PATH_CONN}/adj_matrix_Schaefer_rest_n400_pcor -q1=${PATH_QC}/final_list_subset_qc_task -q2=${PATH_QC}/final_list_subset_qc_rest"
"-m1=rest -m2=movie --path_m1=${PATH_CONN}/adj_matrix_Schaefer_rest_n400_pcor --path_m2=${PATH_CONN}/adj_matrix_Schaefer_movie_n400_pcor -q1=${PATH_QC}/final_list_subset_qc_rest -q2=${PATH_QC}/final_list_subset_qc_movie"
"-m1=movie -m2=rest --path_m1=${PATH_CONN}/adj_matrix_Schaefer_movie_n400_pcor --path_m2=${PATH_CONN}/adj_matrix_Schaefer_rest_n400_pcor -q1=${PATH_QC}/final_list_subset_qc_movie -q2=${PATH_QC}/final_list_subset_qc_rest"
"-m1=task -m2=movie --path_m1=${PATH_CONN}/adj_matrix_Schaefer_task_n400_pcor --path_m2=${PATH_CONN}/adj_matrix_Schaefer_movie_n400_pcor -q1=${PATH_QC}/final_list_subset_qc_task -q2=${PATH_QC}/final_list_subset_qc_movie"
"-m1=movie -m2=task --path_m1=${PATH_CONN}/adj_matrix_Schaefer_movie_n400_pcor --path_m2=${PATH_CONN}/adj_matrix_Schaefer_task_n400_pcor -q1=${PATH_QC}/final_list_subset_qc_movie -q2=${PATH_QC}/final_list_subset_qc_task"
::: "-n=visual" "-n=somatomotor" "-n=dorsal_attention" "-n=salience_ventral_attention" "-n=limbic" "-n=frontoparietal" "-n=default_mode" "-n=Whole_brain" "-n=default_mode_limbic" 
> fpc_schaefer_pcor_args.txt