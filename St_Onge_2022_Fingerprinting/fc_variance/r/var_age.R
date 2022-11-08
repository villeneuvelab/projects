#!/usr/bin/env Rscript
###############################
###############################

#"
## Fingerprinting - CAMCAN Cohort
## Data Analysis - Polynomial (Quadratic) regression - Variability in FC and HCTSA

#Script - Data Analysis - Quadratic regression - Variability
#Author: Frederic St-Onge

#Date created: 2020-12-23
#Date modified: 2021-12-15

#Version: 2.0

#--------
#Purpose:
#--------
#    The purpose of this script is to have a ready script where we can generate results
#    simply and in the same way regardless of the modality used. 
#    
#    In the first iteration, this script generated 1 graph per network, for each modality pair,
#    and it was assembled manually in Power Point. It was a pain. So I remade the script to run from
#    the commandline, but it only need specification on the path to the data, the modalities of interest
#    and the correlation used for the FC matrices. That said, another function is needed if you want 
#    other graphs (e.g., without the polynomial relationship). Could use the "corr" argument in an "if/else" statement.
#   
#
#--------
#Usage:
#--------
#   The script needs to be copied and modifed for each different iteration. What needs to be specifically modified is:
#
#--------
#Background:
#--------
#    (TO WRITE)
#
#    
#--------
#Packages:
#--------
#    The following packages need to be installed:
#      - tidyverse (v. 1.3.0)
#      - glue (v. 1.4.2)
#      - argparser (v. 0.6)
#
#
#--------
#Future updates:
#--------
#    - Find a way to allow for short forms in the parser? Didn't work with the "short=" argument.
#    - Find a way to print stuff to the terminal to indicate help the user determine what is done.
#    - Change the selection process. If the script needs a "network" argument, the pivoting is redundant.
#    - Simplify functions (1 function 1 purpose). Maybe classes?
#
#

################# Step 1 - Laoding libraries and setting the parser

#R libraries (We force messages to be suppressed)

suppressMessages(library(tidyverse)) #General clean-up + ggplot2
suppressMessages(library(glue)) #For glue strings (using arguments)
suppressMessages(library(argparser)) #For CL arguments
suppressMessages(library(patchwork)) #To assemble final graphs
suppressMessages(library(boot)) #For bootstrap
options(scipen = 999) #Replace scientific notation with rounded numbers.

#Argument parser

a <- arg_parser("Generate graphs for Polynomial regression (FPC)")
#cat(a)

a <- add_argument(a, "--path_file", help="Path where the .csv file containing the data to analyse is located.")
a <- add_argument(a, "--output", help="Path where R should export the figure.")
a <- add_argument(a, "--modality", help="Name of the modality 1. Should match columns from file." )
a <- add_argument(a, "--atlas", help="Atlas used to generate the FP. Used to select correct columns and correct output name.")
a <- add_argument(a, "--corr", help="Type of correlation used to generate the matrix. Used to select correct columns and correct output name.")

#Parse the arguments and set to variables

args <- parse_args(a)

path_file <- args$path_file
output <- args$output
modality <- args$modality
atlas <- args$atlas
corr <- args$corr


###################### Step 2 - Prep functions for the analyses

#Data import
import_data <- function(path_file) {
  "Simple function importing the .csv file and returns it."
  
  df <- read_csv(glue("{path_file}"), show_col_types = FALSE)
  #cat("Imported data successfully.")
  
  return(df)
}

#Generate individual graphs
var_plots <- function(data, type, type_label, modality, atlas, network, corr, col_net) {
  "This function plots the FPC graphs and saves it to an object."
  
  p <- data %>% #Imports the data
    ggplot(aes_string(x=glue("age"), y=glue("{type}_{modality}_{atlas}_{network}_{corr}"))) + #Use aes_string to read glue strings
    geom_point(colour = glue("{col_net}"), position="jitter") + #Scatterplot for variance metric
    stat_smooth(method = "lm", formula = y ~ x, colour = glue("{col_net}"), fill = glue("{col_net}")) + #Lm line for variance metric.
    theme_classic() + #Imposes the classic theme (only axis lines)
    theme(panel.grid.major = element_line(), #Reimports the axis grid (1/2)
          panel.grid.minor = element_line(), #Reimports the axis gris (1/2)
          text = element_text(size=14, family="sans"), #Imposes font family and size
          axis.title.x = element_text(face="bold"), #Bolds the axis title
          axis.title.y = element_text(face="bold")) + #Bolds the axis title
    xlab("Age") + #Write labels for X
    ylab(glue("{type_label}")) #Writes labels for Y
  
  #We also want to return a "combined" plot. This is to make it easier to visualize.
  p_comb <- data %>%
    ggplot(aes_string(y=glue("{type}_{modality}_{atlas}_Whole_brain_{corr}"), x=glue("age"))) +
    geom_point(colour = "black", position="jitter", alpha = 0.3) + #Scatterplot for variance metric
    stat_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "black") + #Lm line for variance metric.
    geom_point(data=data, aes_string(y=glue("{type}_{modality}_{atlas}_default_mode_{corr}"), x=glue("age")), colour="#CD3E4E", position = "jitter", alpha = 0.3) +
    stat_smooth(data=data, aes_string(y=glue("{type}_{modality}_{atlas}_default_mode_{corr}"), x=glue("age")), method = "lm", formula = y ~ x, colour = "#CD3E4E", fill = "#CD3E4E") +
    geom_point(data=data, aes_string(y=glue("{type}_{modality}_{atlas}_frontoparietal_{corr}"), x=glue("age")), colour="#E69422", position = "jitter", alpha = 0.3) +
    stat_smooth(data=data, aes_string(y=glue("{type}_{modality}_{atlas}_frontoparietal_{corr}"), x=glue("age")), method = "lm", formula = y ~ x, colour = "#E69422", fill = "#E69422") +
    geom_point(data=data,  aes_string(y=glue("{type}_{modality}_{atlas}_salience_ventral_attention_{corr}"), x=glue("age")), colour="#C43AFA", position = "jitter", alpha = 0.3) +
    stat_smooth(data=data,  aes_string(y=glue("{type}_{modality}_{atlas}_salience_ventral_attention_{corr}"), x=glue("age")), method = "lm", formula = y ~ x, colour = "#C43AFA", fill = "#C43AFA") +
    geom_point(data=data,  aes_string(y=glue("{type}_{modality}_{atlas}_dorsal_attention_{corr}"), x=glue("age")), colour="#00760E", position = "jitter", alpha = 0.3) +
    stat_smooth(data=data,  aes_string(y=glue("{type}_{modality}_{atlas}_dorsal_attention_{corr}"), x=glue("age")), method = "lm", formula = y ~ x, colour = "#00760E", fill = "#00760E") +
    geom_point(data=data,  aes_string(y=glue("{type}_{modality}_{atlas}_limbic_{corr}"), x=glue("age")), colour="#C2B280", position = "jitter", alpha = 0.3) +
    stat_smooth(data=data,  aes_string(y=glue("{type}_{modality}_{atlas}_limbic_{corr}"), x=glue("age")), method = "lm", formula = y ~ x, colour = "#C2B280", fill = "#C2B280") +
    geom_point(data=data,  aes_string(y=glue("{type}_{modality}_{atlas}_somatomotor_{corr}"), x=glue("age")), colour="#4682B4", position = "jitter", alpha = 0.3) +
    stat_smooth(data=data,  aes_string(y=glue("{type}_{modality}_{atlas}_somatomotor_{corr}"), x=glue("age")), method = "lm", formula = y ~ x, colour = "#4682B4", fill = "#4682B4") +
    geom_point(data=data,  aes_string(y=glue("{type}_{modality}_{atlas}_visual_{corr}"), x=glue("age")), colour="#781286", position = "jitter", alpha = 0.3) +
    stat_smooth(data=data,  aes_string(y=glue("{type}_{modality}_{atlas}_visual_{corr}"), x=glue("age")), method = "lm", formula = y ~ x, colour = "#781286", fill = "#781286") +
    theme_classic() + #Imposes the classic theme (only axis lines)
    theme(panel.grid.major = element_line(), #Reimports the axis grid (1/2)
          panel.grid.minor = element_line(), #Reimports the axis gris (1/2)
          text = element_text(size=14, family="sans"), #Imposes font family and size
          axis.title.x = element_text(face="bold"), #Bolds the axis title
          axis.title.y = element_text(face="bold")) + #Bolds the axis title
    xlab("Age") + #Write labels for X
    ylab(glue("{type_label}"))
  
  #Saves the combined plot
  ggsave(glue("{output}/{type}_age_combined_{modality}_{atlas}_{corr}.png"), plot = p_comb, dpi = 700, width = 21, height = 13.5, units = "cm")
  
  
  return(p)
} 

#Generate list of plots based on our task pair modalities
plot_list_gen <- function(data, type, type_label, modality, atlas, corr) {
  "Generates a list of plots, based on a dictionary. Dictionary is based on Yeo networks for now."
  
  #Create dict
  networks_colors <- c(
    "Whole_brain"="black",
    "frontoparietal"="#E69422",
    "salience_ventral_attention"="#C43AFA",
    "default_mode"="#CD3E4E",
    "dorsal_attention"="#00760E",
    "limbic"="#C2B280",
    "somatomotor"="#4682B4",
    "visual"="#781286")
  
  #Create empty list to store the plots
  p <- list()
  
  #Create loop to call the plotting function and save graphs iteratively
  for (i in seq(1, length(networks_colors))){ #For each pair in the dictionary, generate a number. For each number (equal to the number of pairs in the dict), do the following:
    col_net <- networks_colors[i] #Extract the HEX code of the color to use ("value" of the dict)
    network <- names(networks_colors)[i] #Extract the network name ("key" of the dict)
    writeLines(glue("      Drawing plot of network {network}\n"))
    
    p[[i]] <- var_plots(data=data, type=type, type_label=type_label, modality=modality, atlas=atlas, network=network, corr=corr, col_net=col_net)}
  
  return(p) #Returns a list of plots
}

#Generate the patchwork
patch_sew <- function(plots_list, type, output, modality, atlas, corr){
  "Simple function generating a patchwork. It takes the output of 'plot_list_gen' as its input. It also saves the "
  
  sewn_patch <- wrap_plots(plots_list) #It should arrange the plots in 2 rows of 3 and 1 row of 2.
  
  ggsave(glue("{output}/{type}_{modality}_{atlas}_{corr}.png"), plot = sewn_patch, dpi = 700, width = 21, height = 13.5, units = "in")
  #Saves the figure to a file, based on the strings used to generate the initial figure. Haven't figured out the best dimensions for the overall plot yet.
  
}

#Prep bootstrap functions for the linear models
boot_coefs <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data=d) #Uses the formula with set 
  return(coef(fit))
}
boot_rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data=d)
  return(summary(fit)$adj.r.square)
}

#Linear models calculation and export
lm_func <- function(data, type, modality, atlas, network, corr) {
  "Function where we compute different modelling (i.e., with age only, with covariates (backward elim), with bootstrap on both"
  set.seed(667)
  writeLines(glue("    Working on {type}, {modality}, {atlas}, {network}, {corr}\n"))
  
  #Simple linear model. Need 0 NA for bootstrap.
  #data_clean_simple <- data %>%
  #  select(glue("{type}_{modality}_{atlas}_{network}_{corr}"), age) %>%
  #  na.omit()
  
  #lm_reg <- lm(glue("{type}_{modality}_{atlas}_{network}_{corr} ~ age"), data=data_clean_simple)
  
  #Bootstrap the coefficient and the adjR2
  #boot_lm_coefs <- boot(data=data, statistic = boot_coefs, R=1000, formula=glue("{type}_{modality}_{atlas}_{network}_{corr} ~ age"))
  #boot_lm_rsq <- boot(data=data, statistic = boot_rsq, R=1000, formula=glue("{type}_{modality}_{atlas}_{network}_{corr} ~ age"))
  
  #Backward elimination model. Again, need 0 NA for bootstrap.
  data_clean_full <- data %>%
    select(glue("{type}_{modality}_{atlas}_{network}_{corr}"), age, handedness, sex_string,
           glue("fd_scr_{modality}"), glue("nb_frames_{modality}")) %>%
    na.omit()
  
  #Full model with covariates (we just prep the formula here)
  full_formula <- glue("{type}_{modality}_{atlas}_{network}_{corr} ~ age + handedness + sex_string + 
                      fd_scr_{modality} + nb_frames_{modality}")
  
  lm_full <- lm(full_formula, data=data_clean_full)
  
  #Bootstrap the coefficient and adj R2
  boot_full_coefs <- boot(data=data, statistic = boot_coefs, R=1000, formula=full_formula)
  boot_full_rsq <- boot(data=data, statistic = boot_rsq, R=1000, formula=full_formula)
  
  #Prep the prints. We want the full models followed
  sink(glue("{output}/lm_boot_{type}_{network}_{modality}_{atlas}_{corr}.txt"))
  print("###################")
  print(glue("----- Quadratic regression: {type} ~ age ({network}, {modality}, {atlas}, {corr})"))
  #print("--- Simple linear regression")
  #print(summary(lm_reg))
  #print(" ")
  #print("--------------------------------")
  #print("--- Bootstrapped coefficient (age and age^2) - Simple linear regression")
  #print(boot.ci(boot_lm_coefs, type="bca", index=2)$bca)
  #print(" ")
  #print("--------------------------------")
  #print("--- Bootstrapped determination coefficient (R2) - Simple linear regression")
  #print(boot.ci(boot_lm_rsq, type="bca")$bca)
  #print(" ")
  print("--------------------------------")
  print("--- Full model with covariates")
  print(summary(lm_full))
  print(" ")
  print("--------------------------------")
  print("--- Bootstrapped coefficent (age and age^2 - Multiple regression model")
  print(boot.ci(boot_full_coefs, type="bca", index=2)$bca)
  print(" ")
  print("--------------------------------")
  print("--- Bootstrapped determination coefficient (R2) - Simple linear regression")
  print(boot.ci(boot_full_rsq, type="bca")$bca)
  print("###################")
  sink()
}

#Loop to run the LMs for each networks
loop_lm <- function(data, type, modality, atlas, corr) {
  
  #Create list
  net_list <- c("Whole_brain", "frontoparietal", "salience_ventral_attention",
                "default_mode", "dorsal_attention", "limbic", "somatomotor", "visual")
  
  for (net in net_list) {
    
    lm_func(data=data, type=type, modality=modality, atlas=atlas, network=net, corr=corr)
    
  }
  
}

###################### Step 3 - Launch functions

writeLines(glue("###################################\n"))
writeLines(glue("Launching Var ~ Age script analysis\n"))
writeLines(glue("   - Parameters:\n"))
writeLines(glue("       - .csv file used: {path_file}\n"))
writeLines(glue("       - Output will be here: {output}\n"))
writeLines(glue("\n"))
writeLines(glue("       - Modality: {modality}\n"))
writeLines(glue("       - Atlas used: {atlas}\n"))
writeLines(glue("       - Correlation: {corr}\n"))
writeLines(glue("-----------------------------------\n"))
writeLines(glue("Importing the data...\n"))
data <- import_data(path_file = glue("{path_file}"))
writeLines(glue("Generating the plots...\n"))
p1 <- plot_list_gen(data=data, type="var", type_label="Variance coefficient", modality=glue("{modality}"), atlas=glue("{atlas}"), corr=glue({"{corr}"}))
p2 <- plot_list_gen(data=data, type="mean_dist", type_label="Mean absolute distance", modality=glue("{modality}"), atlas=glue("{atlas}"), corr=glue({"{corr}"}))
writeLines(glue("Sewing the plots together and saving to file...\n"))
patch_sew(plots_list = p1, type="var", output = glue("{output}"), modality = glue("{modality}"), atlas=glue("{atlas}"), corr=glue("{corr}"))
patch_sew(plots_list = p2, type="mean_dist", output = glue("{output}"), modality = glue("{modality}"), atlas=glue("{atlas}"), corr=glue("{corr}"))

writeLines(glue("Launching the linear models for var_coef...\n"))
loop_lm(data=data, type="var", modality=glue("{modality}"), atlas=glue("{atlas}"), corr=glue("{corr}"))
writeLines(glue("Launching the linear models for mean_dist_var...\n"))
loop_lm(data=data, type="mean_dist", modality=glue("{modality}"), atlas=glue("{atlas}"), corr=glue("{corr}"))
writeLines("Done!")



