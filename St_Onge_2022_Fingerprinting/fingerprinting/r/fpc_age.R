#!/usr/bin/env Rscript
###############################
###############################

#"
## Fingerprinting - CAMCAN Cohort
## Data Analysis - Polynomial (Quadratic) regression

#Script - Data Analysis - Quadratic regression - FPC/ASC
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
suppressMessages(library(lmtest)) #For likelihood ratio test.
options(scipen = 999) #Replace scientific notation with rounded numbers.

#Argument parser

a <- arg_parser("Generate graphs for Polynomial regression (FPC)")
#cat(a)

a <- add_argument(a, "--path_file", help="Path where the .csv file containing the data to analyse is located.")
a <- add_argument(a, "--output", help="Path where R should export the figure.")
a <- add_argument(a, "--modality1", help="Name of the modality 1. Should match columns from file." )
a <- add_argument(a, "--modality2", help="Name of the modality 2. Should match columns from file." )
a <- add_argument(a, "--atlas", help="Atlas used to generate the FP. Used to select correct columns and correct output name.")
a <- add_argument(a, "--corr", help="Type of correlation used to generate the matrix. Used to select correct columns and correct output name.")
a <- add_argument(a, "--col_asc", default="grey", help="Color used for ASC")
a <- add_argument(a, "--type", help='Type of network used for fingerprinting')

#Parse the arguments and set to variables

args <- parse_args(a)

path_file <- args$path_file
output <- args$output
modality1 <- args$modality1
modality2 <- args$modality2
atlas <- args$atlas
corr <- args$corr
col_asc <- args$col_asc
type <- args$type

###################### Step 2 - Prep functions for the analyses

#Data import
import_data <- function(path_file) {
  "Simple function importing the .csv file and returns it."
  
  df <- read_csv(glue("{path_file}"), show_col_types = FALSE)
  #cat("Imported data successfully.")
  
  return(df)
}

round_p_funct <- function(pval) {
  
  if (pval < 0.001) {
    rounded_p = "< 0.001"
  } else {
    rounded_p = format(round(pval, 2), nsmall = 2)
  }
  
  return(rounded_p)
}

pstar_funct <- function(pval) {
  
  if (pval < 0.001) {
    pstar = "***"
  } else if (pval < 0.01) {
    pstar = "**"
  } else if (pval < 0.05) {
    pstar = "*"
  } else if (pval < 0.1) {
    pstar = "\U00B0"
  } else {
    pstar = ""
  }
  
  return(pstar)
}

#Generate individual graphs
fpc_asc_plots <- function(data, modality1, modality2, atlas, network, corr, col_net, col_asc, type) {
  "This function plots the FPC graphs and saves it to an object."
  
  #Compute the Stimson equation parameters
  lm_simp_fpc <- lm(glue("fpc_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr} ~ age"), data=data)
  lm_quad_fpc <- lm(glue("fpc_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr} ~ age + I(age^2)"), data=data)
  b0_fpc <- summary(lm_quad_fpc)$coefficients[1] #Intercept
  b1_fpc <- summary(lm_quad_fpc)$coefficients[2] #Coefficient of age alone
  b1_fpc_pval <- summary(lm_quad_fpc)$coefficients[2,4]
  b1_fpc_pstars <- pstar_funct(b1_fpc_pval)
  b2_fpc <- summary(lm_quad_fpc)$coefficients[3] #Coefficient of quadratic of age
  b2_fpc_pval <- summary(lm_quad_fpc)$coefficients[3,4]
  b2_fpc_pstars <- pstar_funct(b2_fpc_pval)
  
  m_fpc <- (b0_fpc - ((b1_fpc^2) / (4*b2_fpc))) #Minimum/maximum value of Y at inflexion
  f_fpc <- as.integer((-b1_fpc) / (2*b2_fpc)) #Minimum/maximum value of X at inflexion
  
  lm_simp_asc <- lm(glue("asc_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr} ~ age"), data=data)
  lm_quad_asc <- lm(glue("asc_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr} ~ age + I(age^2)"), data=data)
  b0_asc <- summary(lm_quad_asc)$coefficients[1] #Intercept
  b1_asc <- summary(lm_quad_asc)$coefficients[2] #Coefficient of age alone
  b1_asc_pval <- summary(lm_quad_asc)$coefficients[2,4]
  b1_asc_pstars <- pstar_funct(b1_asc_pval)
  b2_asc <- summary(lm_quad_asc)$coefficients[3] #Coefficient of quadratic of age
  b2_asc_pval <- summary(lm_quad_asc)$coefficients[3,4]
  b2_asc_pstars <- pstar_funct(b2_asc_pval)
  
  r2adj_fpc <- format(round(summary(lm_quad_fpc)$adj.r.squared, 2), nsmall=2)
  r2adj_asc <- format(round(summary(lm_quad_asc)$adj.r.squared, 2), nsmall=2)
  
  m_asc <- (b0_asc - ((b1_asc^2) / (4*b2_asc))) #Minimum/maximum value of Y at inflexion
  f_asc <- as.integer((-b1_asc) / (2*b2_asc))
  
  #Likelihood ratio test for Age2
  lr_test_fpc <- lrtest(lm_simp_fpc, lm_quad_fpc)
  pval_lr_test_fpc <- lr_test_fpc$`Pr(>Chisq)`[[2]]
  r_pval_lr_test_fpc <- round_p_funct(pval_lr_test_fpc)
  
  lr_test_asc <- lrtest(lm_simp_asc, lm_quad_asc)
  pval_lr_test_asc <- lr_test_asc$`Pr(>Chisq)`[[2]]
  r_pval_lr_test_asc <- round_p_funct(pval_lr_test_asc)
  
  #Rounded betas
  b1_fpc_round <- format(b1_fpc, scientific=TRUE, digits=3)
  b2_fpc_round <- format(b2_fpc, scientific=TRUE, digits=3)
  b1_asc_round <- format(b1_asc, scientific=TRUE, digits=3)
  b2_asc_round <- format(b2_asc, scientific=TRUE, digits=3)
  
  
  #Plot the data
  p <- data %>% #Imports the data
    ggplot(aes_string(x=glue("age"), y=glue("fpc_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr}"))) + #Use aes_string to read glue strings
    geom_point(colour = glue("{col_net}"), position="jitter", size=1) + #Scatterplot for FPC
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour = glue("{col_net}"), fill = glue("{col_net}")) + #Lm line for FPC
    {if(f_fpc > 1 & f_fpc < 100)annotate('segment', x = f_fpc, xend = f_fpc, y = m_fpc*1.03, yend = m_fpc + 0.13, size = 0.4, color = glue("{col_net}"))} + #If the age value "makes sense" (i.e., not extreme), we can plot the inflexion point
    {if(f_fpc > 1 & f_fpc < 100)annotate('text', label=glue("{f_fpc} years"), x=f_fpc + 5.5, y=m_fpc+0.125, color = "black", alpha = 0.8, size=3)} + #If the age value "makes sense" (i.e., not extreme), we can plot the inflexion point.
    geom_point(data = data, aes_string(x=glue("age"), y=glue("asc_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr}")), colour = glue("{col_asc}"), position = "jitter", size=1) + #Scatterplot for ASC (needs to respecify the data and aesthetic)
    stat_smooth(data = data, aes_string(x=glue("age"), y=glue("asc_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr}")), method = "lm", formula = y ~ x + I(x^2), colour = glue("{col_asc}"), fill = glue("{col_asc}")) + #Lm line for ASC (needs to respecify the data and aesthetic)
    {if(f_asc > 1 & f_asc < 100)annotate('segment', x=f_asc, xend=f_asc, y=m_asc*0.97, yend = m_asc - 0.05, size=0.4, color=glue("{col_asc}"))} +
    {if(f_asc > 1 & f_asc < 100)annotate('text',label=glue("{f_asc} years"), x=f_asc + 5.5, y=m_asc-0.045, color="black", alpha = 0.8, size=3)} +
    annotate("text", y=Inf, x =Inf, hjust=1, vjust=1,  
             label=glue("\U03B2 age = {b1_fpc_round} {b1_fpc_pstars}\n\U03B2 age\U00B2 = {b2_fpc_round} {b2_fpc_pstars}\nR\U00B2 adj = {r2adj_fpc}\nLR Test = {r_pval_lr_test_fpc}"),
                        size=3) +
    annotate("text", y=-Inf, x =-Inf, hjust=-0.1, vjust=-0.1,  
             label=glue("\U03B2 age = {b1_asc_round} {b1_asc_pstars}\n\U03B2 age\U00B2 = {b2_asc_round} {b2_asc_pstars}\nR\U00B2 adj = {r2adj_asc}\nLR Test = {r_pval_lr_test_asc}"),
                        size=3) +
    theme_classic() + #Imposes the classic theme (only axis lines)
    theme(panel.grid.major = element_line(), #Reimports the axis grid (1/2)
          panel.grid.minor = element_line(), #Reimports the axis gris (1/2)
          text = element_text(size=12, family="sans"), #Imposes font family and size
          axis.title.x = element_text(face="bold"), #Bolds the axis title
          axis.title.y = element_text(face="bold")) + #Bolds the axis title
    xlab("Age") + #Write labels for X
    ylab("Correlation coefficient") #Writes labels for Y
  
  return(p)
}

#Generate list of plots based on our task pair modalities
plot_list_gen <- function(data, modality1, modality2, atlas, corr, type) {
  "Generates a list of plots, based on a dictionary. Dictionary is based on Yeo networks for now."
  
  #Create dict
  if (type == "within" & atlas == "Schaefer") {
    networks_colors <- c(
      "Whole_brain"="black",
      "frontoparietal"="#E69422",
      "salience_ventral_attention"="#C43AFA",
      "default_mode"="#CD3E4E",
      "dorsal_attention"="#00760E",
      "limbic"="#C2B280",
      "somatomotor"="#4682B4",
      "visual"="#781286",
      "random1"="#f4c2c2",
      "random2"="brown")
  } else if (type == "between" & atlas == 'Schaefer') {
    networks_colors <- c(
      "frontoparietal"="#E69422",
      "salience_ventral_attention"="#C43AFA",
      "default_mode"="#CD3E4E",
      "dorsal_attention"="#00760E",
      "limbic"="#C2B280",
      "somatomotor"="#4682B4",
      "visual"="#781286")
  } else if (type == "within" & atlas == "Power") {
      networks_colors <- c(
          "uncertain"="#FFF2CC", #Beige
          "salience"="#5D3FD3", #Dark purple (instead of black)
          "cingulo_opercular"="#800080", #Purple
          "dorsal_attention"="#008000", #Green
          "frontoparietal"="#CCCC00", #Dark Yellow 
          "auditory"="#FFC0CB", #Pink
          "default_mode"="#FF0000", #Red
          "visual"="#0000FF", #Blue
          "cerebellar"="#9FC5E8", #Pale blue
          "ventral_attention"="#008080", #Cyan
          "sensory_somatomotor"="#00FFFF", #Orange
          "Whole_brain"="black",
          "subcortical"="#A52A2A") #Brown
    }

  
  #Create empty list to store the plots
  p <- list()
  
  #Create loop to call the plotting function and save graphs iteratively
  for (i in seq(1, length(networks_colors))){ #For each pair in the dictionary, generate a number. For each number (equal to the number of pairs in the dict), do the following:
    col_net <- networks_colors[i] #Extract the HEX code of the color to use ("value" of the dict)
    network <- names(networks_colors)[i] #Extract the network name ("key" of the dict)
    writeLines(glue("      Drawing plot of network {network}"))
    
    p[[i]] <- fpc_asc_plots(data=data, modality1=modality1, modality2=modality2, atlas=atlas, network=network, corr=corr, col_net=col_net, col_asc=col_asc, type=type)}
  
  return(p) #Returns a list of plots
}

#Generate the patchwork
patch_sew <- function(plots_list, output, modality1, modality2, atlas, corr, type){
  "Simple function generating a patchwork. It takes the output of 'plot_list_gen' as its input. It also saves the "
  
  sewn_patch <- wrap_plots(plots_list, ncol=3) #It should arrange the plots in 2 rows of 3 and 1 row of 2.
  
  ggsave(glue("{output}/fpc_asc_{type}_{modality1}_{modality2}_{atlas}_{corr}.png"), plot = sewn_patch, dpi = 700, width = 35.6, height=35.6, units = "cm")
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
lm_func <- function(data, metric, modality1, modality2, atlas, network, corr, type) {
  "Function where we compute different modelling (i.e., with age only, with covariates (backward elim), with bootstrap on both"
  set.seed(667)
  writeLines(glue("    Working on {metric}, {modality1}, {modality2}, {type}, {atlas}, {network}, {corr}\n"))
  
  #Simple linear model. Need 0 NA for bootstrap.
  #data_clean_simple <- data %>%
  #  select(glue("{metric}_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr}"), age) %>%
  #  na.omit()
  
  #lm_reg <- lm(glue("{metric}_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr} ~ age + I(age^2)"), data=data_clean_simple)
  
  #Bootstrap the coefficient and the adjR2
  #boot_lm_coefs <- boot(data=data, statistic = boot_coefs, R=1000, formula=glue("{metric}_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr} ~ age + I(age^2)"))
  #boot_lm_rsq <- boot(data=data, statistic = boot_rsq, R=1000, formula=glue("{metric}_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr} ~ age + I(age^2)"))
  
  #Backward elimination model. Again, need 0 NA for bootstrap.
  data_clean_full <- data %>%
    select(glue("{metric}_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr}"), age, handedness, sex_string,
           glue("fd_scr_{modality1}"), glue("fd_scr_{modality2}"), glue("nb_frames_{modality1}"), 
           glue("nb_frames_{modality2}")) %>%
    na.omit()
  
  #Full model with covariates (we just prep the formula here)
  full_formula <- glue("{metric}_{modality1}_{modality2}_{type}_{atlas}_{network}_{corr} ~ age + I(age^2) + handedness + sex_string + 
                      fd_scr_{modality1} + fd_scr_{modality2} + nb_frames_{modality1} + nb_frames_{modality2}")
  
  lm_full <- lm(full_formula, data=data_clean_full)
  
  #Bootstrap the coefficient and adj R2
  boot_full_coefs <- boot(data=data, statistic = boot_coefs, R=1000, formula=full_formula)
  boot_full_rsq <- boot(data=data, statistic = boot_rsq, R=1000, formula=full_formula)
  
  #Prep the prints. We want the full models followed
  sink(glue("{output}/lm_boot_{metric}_{type}_{network}_{modality1}_{modality2}_{atlas}_{corr}.txt"))
  print("###################")
  print(glue("----- Quadratic regression: {metric} ~ age ({type}, {network}, {modality1}, {modality2}, {atlas}, {corr})"))
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
loop_lm <- function(data, metric, modality1, modality2, atlas, corr, type) {
  
  #Create list
  net_list = NULL
  if (type == "within" & atlas == "Schaefer"){
    net_list <- c("Whole_brain", "frontoparietal", "salience_ventral_attention",
                  "default_mode", "dorsal_attention", "limbic", "somatomotor", "visual", "random1", "random2")
  } else {
    net_list <- c("frontoparietal", "salience_ventral_attention",
                  "default_mode", "dorsal_attention", "limbic", "somatomotor", "visual")
  }
  
  for (net in net_list) {
        lm_func(data=data, metric=metric, modality1=modality1, modality2=modality2, atlas=atlas, network=net, corr=corr, type=type)
      }
  
}


###################### Step 3 - Launch functions

writeLines(glue("###################################\n"))
writeLines(glue("Launching Age ~ FPC script analysis\n"))
writeLines(glue("   - Parameters:\n"))
writeLines(glue("       - .csv file used: {path_file}\n"))
writeLines(glue("       - Output will be here: {output}\n"))
writeLines(glue("\n"))
writeLines(glue("       - Modality1: {modality1}\n"))
writeLines(glue("       - Modality2: {modality2}\n"))
writeLines(glue("       - Atlas used: {atlas}\n"))
writeLines(glue("       - Correlation: {corr}\n"))
writeLines(glue("       - Type: {type}\n"))
writeLines(glue("-----------------------------------\n"))
writeLines(glue("Importing the data...\n"))
data <- import_data(path_file = glue("{path_file}"))
writeLines(glue("Generating the plots..."))
p <- plot_list_gen(data=data, modality1=glue("{modality1}"), modality2=glue("{modality2}"), atlas=glue("{atlas}"), corr=glue({"{corr}"}), type=glue({"{type}"}))
writeLines(glue("Sewing the plots together and saving to file..."))
patch_sew(plots_list = p, output = glue("{output}"), modality1 = glue("{modality1}"), modality2 = glue("{modality2}"), atlas=glue("{atlas}"), corr=glue("{corr}"), type=glue("{type}"))

writeLines(glue("Launching the linear models for FPC..."))
loop_lm(data=data, metric="fpc", modality1=glue("{modality1}"), modality2=glue("{modality2}"), atlas=glue("{atlas}"), corr=glue("{corr}"), type=glue("{type}"))
writeLines(glue("Launching the linear models for ASC..."))
loop_lm(data=data, metric="asc", modality1=glue("{modality1}"), modality2=glue("{modality2}"), atlas=glue("{atlas}"), corr=glue("{corr}"), type=glue("{type}"))
writeLines("Done!")


















