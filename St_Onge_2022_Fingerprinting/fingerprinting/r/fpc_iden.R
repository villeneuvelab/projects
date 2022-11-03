#!/usr/bin/env Rscript

#"
## Fingerprinting - CAMCAN Cohort
## Data Analysis - Objective 1 - Describing the fingerprinting across the lifespan

#Script - Data Analysis - Objective 1
#Author: Frédéric St-Onge

#Date created: 2020-12-17
#Date modified: 2022-01-26

#Updates:
### 2022-01-26: Added compatibility to generate graphs for "between" network connectivity.

#Version: 1.1

#--------
#Purpose:
#--------
#    The purpose of this script is to have a command-line ready script to generate the analyses quickly
#    once tested in a Markdown document.
#    
#    In this first objective, we want to:
#      1. Generate tables where the percentage of Fingerprint success/failures is reported with confidence intervals,
#      for each network within each pair of subjects.
#      2. Generate bar plots to illustrate these percentages.
#   
#    Note: The script is still quite rigid, and will only work for the Schaefer atlas. This is because the functions
#    generate graphs specific to the networks of interest. This will need to be changed when we will do the 
#    Power atlas. 
#    
#    It will also calculate the pairs of modalities based on what we have in Cam-CAN
#
#--------
#Usage:
#--------
#   The script needs to be copied and modifed for each different iteration. What needs to be specifically modified is:
#   - The file used to generate the stats (file that was cleaned in the Data Cleaning step)
#   - The 'atlas' used to generate the fingerprinting (but now only works in Schaefer)
#   - The 'correlation' used to generate the fingerprinting
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
#
#"

########## Step 1 - Setting up the script

#Loading R packages
suppressMessages(library(tidyverse)) #To tidy the dataframes before generating the number and plot
suppressMessages(library(glue)) #To attach arguments to functions we run and to generate custom names
library(argparser) #To feed arguments from the command line

#Setting the arguments

a <- arg_parser("Generates Table and Graph for Objective 1 (1 pair of modality)")

a <- add_argument(a, "--path_file", help="Path where the .csv file containing the data to analyse is located.")
a <- add_argument(a, "--output", help="Path where R should export the table and figure.")
a <- add_argument(a, "--modality1", help="Name of the modality 1. Should match columns from file." )
a <- add_argument(a, "--modality2", help="Name of the modality 2. Should match columns from file." )
a <- add_argument(a, "--atlas", help="Atlas used to generate the FP. Used to select correct columns and correct output name.")
a <- add_argument(a, "--corr", help="Type of correlation used to generate the matrix. Used to select correct columns and correct output name.")
a <- add_argument(a, "--type", help="Whether the FPC in use is within or between-network")

#Parse the arguments and set to variables
args <- parse_args(a)

path_file <- args$path_file
output <- args$output
modality1 <- args$modality1
modality2 <- args$modality2
atlas <- args$atlas
corr <- args$corr
type <- args$type
#writeLines(glue("Showing values for: {path_file} {output} {modality1} {modality2} {atlas} {corr} {type}"))
#Importing the .csv file

fpc_data_clean <- read.csv(path_file)

#Setting functions to be used in the script

############
#----------------
# Function to generate the percentage calculation
#----------------
percentage_fingerprinted_calculator <- function(df, modality1, modality2, atlas, corr, type) {
  
  # Selects the columns of the data that check fingerprinting status for the correct modality
  perc_df <- NULL #Creates an empty variable that will be replaced depending on the if/else
  
  if (type == "within" & atlas == "Schaefer") {
    perc_df <- df %>%
      select(id_camcan, contains(glue('non_fp_{modality1}_{modality2}_{type}'))) %>%
      na.omit() %>%
      rename(
        whole_brain = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_Whole_brain_{corr}"),
        frontopar = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_frontoparietal_{corr}"),
        sal_vent = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_salience_ventral_attention_{corr}"),
        dmn = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_default_mode_{corr}"),
        dor_att = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_dorsal_attention_{corr}"),
        lim = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_limbic_{corr}"),
        somat = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_somatomotor_{corr}"),
        vis = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_visual_{corr}"),
        rand1 = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_random1_{corr}"),
        rand2 = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_random2_{corr}"))
    
  } else if (type == "between" & atlas == "Schaefer") {
    perc_df <- df %>%
      select(id_camcan, contains(glue('non_fp_{modality1}_{modality2}_{type}'))) %>%
      na.omit() %>%
      rename(
        frontopar = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_frontoparietal_{corr}"),
        sal_vent = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_salience_ventral_attention_{corr}"),
        dmn = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_default_mode_{corr}"),
        dor_att = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_dorsal_attention_{corr}"),
        lim = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_limbic_{corr}"),
        somat = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_somatomotor_{corr}"),
        vis = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_visual_{corr}"))
    
  } else if (type == "within" & atlas == "Power") {
      perc_df <- df %>%
          select(id_camcan, contains(glue('non_fp_{modality1}_{modality2}_{type}'))) %>%
          na.omit() %>%
          rename(
              uncert = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_uncertain_{corr}"),
              salience = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_salience_{corr}"),
              cing_opp = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_cingulo_opercular_{corr}"),
              dors_att = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_dorsal_attention_{corr}"),
              frontopar = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_frontoparietal_{corr}"),
              audit = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_auditory_{corr}"),
              dmn = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_default_mode_{corr}"),
              vis = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_visual_{corr}"),
              cereb = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_cerebellar_{corr}"),
              vent = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_ventral_attention_{corr}"),
              somat = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_sensory_somatomotor_{corr}"),
              whole_brain = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_Whole_brain_{corr}"),
              sub = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_subcortical_{corr}")
          )
  }

  
  head(perc_df)
  # Creates a dataframe where the we calculate the percentage for each network
  
  perc_calc <- NULL
  
  if (type == "within" & atlas == "Schaefer") {
    perc_calc <- data.frame(
      network = c("Whole Brain", "Frontoparietal", "Salience/ Ventral attention", "Default-Mode", "Dorsal Attention", "Limbic", 
                  "Somatomotor","Visual", "Random1", "Random2"),
      
      percentage_recon = c((sum(perc_df$whole_brain)/length(perc_df$whole_brain)), 
                           (sum(perc_df$frontopar)/length(perc_df$frontopar)), 
                           (sum(perc_df$sal_vent)/length(perc_df$sal_vent)),
                           (sum(perc_df$dmn)/length(perc_df$dmn)),
                           (sum(perc_df$dor_att)/length(perc_df$dor_att)),
                           (sum(perc_df$lim)/length(perc_df$lim)),
                           (sum(perc_df$somat)/length(perc_df$somat)),
                           (sum(perc_df$vis)/length(perc_df$vis)),
                           (sum(perc_df$rand1)/length(perc_df$rand1)),
                           (sum(perc_df$rand2)/length(perc_df$rand2))
                           ),
      
      number_pos = c(sum(perc_df$whole_brain), 
                     sum(perc_df$frontopar), 
                     sum(perc_df$sal_vent),
                     sum(perc_df$dmn),
                     sum(perc_df$dor_att),
                     sum(perc_df$lim),
                     sum(perc_df$somat),
                     sum(perc_df$vis),
                     sum(perc_df$rand1),
                     sum(perc_df$rand2)
                     ),
      
      total = c(length(perc_df$whole_brain), 
                length(perc_df$frontopar), 
                length(perc_df$sal_vent),
                length(perc_df$dmn),
                length(perc_df$dor_att),
                length(perc_df$lim),
                length(perc_df$somat),
                length(perc_df$vis),
                length(perc_df$rand1),
                length(perc_df$rand2))
      
    )
  } else if (type == "between" & atlas == "Schaefer") {
    perc_calc <- data.frame(
      network = c("Frontoparietal", "Salience/ Ventral attention", "Default-Mode", "Dorsal Attention", "Limbic", "Somatomotor","Visual"),
      
      percentage_recon = c((sum(perc_df$frontopar)/length(perc_df$frontopar)), 
                           (sum(perc_df$sal_vent)/length(perc_df$sal_vent)),
                           (sum(perc_df$dmn)/length(perc_df$dmn)),
                           (sum(perc_df$dor_att)/length(perc_df$dor_att)),
                           (sum(perc_df$lim)/length(perc_df$lim)),
                           (sum(perc_df$somat)/length(perc_df$somat)),
                           (sum(perc_df$vis)/length(perc_df$vis))),
      
      number_pos = c(sum(perc_df$frontopar), 
                     sum(perc_df$sal_vent),
                     sum(perc_df$dmn),
                     sum(perc_df$dor_att),
                     sum(perc_df$lim),
                     sum(perc_df$somat),
                     sum(perc_df$vis)),
      
      total = c(length(perc_df$frontopar), 
                length(perc_df$sal_vent),
                length(perc_df$dmn),
                length(perc_df$dor_att),
                length(perc_df$lim),
                length(perc_df$somat),
                length(perc_df$vis))
      
    )
  } else if (type == "within" & atlas == "Power") {
      perc_calc <- data.frame(
          network = c("Uncertain", "Salience", "Cingular-Oppercular", "Dorsal Attention", "Frontoparietal", "Auditory",
                      "Default-mode", "Visual", "Cerebellar", "Ventral attention", "Sensory Somatomotor", "Whole brain", "Subcortical"),
          percentage_recon = c((sum(perc_df$uncert)/length(perc_df$uncert)), 
                               (sum(perc_df$salience)/length(perc_df$salience)),
                               (sum(perc_df$cing_opp)/length(perc_df$cing_opp)),
                               (sum(perc_df$dors_att)/length(perc_df$dors_att)),
                               (sum(perc_df$frontopar)/length(perc_df$frontopar)),
                               (sum(perc_df$audit)/length(perc_df$audit)),
                               (sum(perc_df$dmn)/length(perc_df$dmn)),
                               (sum(perc_df$vis)/length(perc_df$vis)),
                               (sum(perc_df$cereb)/length(perc_df$cereb)),
                               (sum(perc_df$vent)/length(perc_df$vent)),
                               (sum(perc_df$somat)/length(perc_df$somat)),
                               (sum(perc_df$whole_brain)/length(perc_df$whole_brain)),
                               (sum(perc_df$sub)/length(perc_df$sub))
                               ),
          
          number_pos = c(sum(perc_df$uncert), 
                         sum(perc_df$salience),
                         sum(perc_df$cing_opp),
                         sum(perc_df$dors_att),
                         sum(perc_df$frontopar),
                         sum(perc_df$audit),
                         sum(perc_df$dmn),
                         sum(perc_df$vis),
                         sum(perc_df$cereb),
                         sum(perc_df$vent),
                         sum(perc_df$somat),
                         sum(perc_df$whole_brain),
                         sum(perc_df$sub)),
          
          total = c(length(perc_df$uncert), 
                    length(perc_df$salience),
                    length(perc_df$cing_opp),
                    length(perc_df$dors_att),
                    length(perc_df$frontopar),
                    length(perc_df$audit),
                    length(perc_df$dmn),
                    length(perc_df$vis),
                    length(perc_df$cereb),
                    length(perc_df$vent),
                    length(perc_df$somat),
                    length(perc_df$whole_brain),
                    length(perc_df$sub)
          )
      )
  }

  #We calculate confidence intervals for each network based on an alpha threshold of 0.05.
  perc_calc$low_conf_interval <- perc_calc$percentage_recon - 1.96*sqrt(perc_calc$percentage_recon*(1-perc_calc$percentage_recon)/perc_calc$total)
  perc_calc$high_conf_interval <- perc_calc$percentage_recon + 1.96*sqrt(perc_calc$percentage_recon*(1-perc_calc$percentage_recon)/perc_calc$total)
  
  return(perc_calc)
  
}
#---------------
# Function to generate the bar graphs
bar_graph_generator <- function(df, modality1, modality2, atlas, corr, type, fill_pass="#024F6D", fill_fail="#ED1B2F") {
  "
  General description:
  --------------------
  This function is made to generate bar graphs for the fingerprinting. It needs as input:
    1) The two modalities of interest, 2) the atlas used and 3) the correlation used for
    the matrices.
    
  As optional input, the color of the fill for the graph can be modified (fill_pass and
  fill_fail options). By default, it comes out as navy blue for the pass and crimson red
  for the fail. Any other changes to the graph requires a modification of the function.
  
  As a general note, the graphing functions used in this script use minimal options. This
  is to make it easier to modify the graph
  outside of R.
    
  General function:
  -----------------
  The script, in order:
    1. Imports the data from the initial dataframe
    2. Selects the non_fp columns of the correct pair of modality and rename columns
    3. Pivots the data from wide to long
    4. Create a new column for whether or not they passed Fingerprinting
    5. Eliminate missing values
    
    6. Plots the stacked bar graph (x = networks, y = number of participants passing 
    (or not) fingerprinting).
  
  "
  
  #Define dataframe for the graph (without missing values)
  p <- NULL
  if (type == "within" & atlas == "Schaefer") {
    df1 <- df %>%
      select(id_camcan, contains(glue('non_fp_{modality1}_{modality2}_{type}'))) %>% #Selecting only the columns with "modality1/modality2"
      rename(
        "Whole Brain" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_Whole_brain_{corr}"),
        "Frontoparietal" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_frontoparietal_{corr}"),
        "Salience Ventral Attention" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_salience_ventral_attention_{corr}"),
        "Default-Mode" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_default_mode_{corr}"),
        "Dorsal Attention" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_dorsal_attention_{corr}"),
        "Limbic" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_limbic_{corr}"),
        "Somatomotor" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_somatomotor_{corr}"),
        "Visual" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_visual_{corr}"),
        "Random 1" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_random1_{corr}"),
        "Random 2" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_random2_{corr}")
      ) %>% #Rename all columns to a easy acronym to avoid messing with levels instead
      pivot_longer(!id_camcan, names_to="network", values_to="identification") %>% #Pivoting the dataframe to a long format (ggplot won't work otherwise)
      mutate(label_pos = factor(case_when(identification == 1 ~ "Identified", identification == 0 ~ "Not identified", TRUE ~ NA_character_))) %>% #Adding a "factor" column for "identified/not_identified" for the "fill" of the the graph part ()
      na.omit() #Omitting rows with missing data
    
    #Define the graph to produce
    p <- df1 %>%
      ggplot(aes(y=network, fill=factor(label_pos, levels=c("Not identified", "Identified")))) + #Bar graph where 1 column = 1 network
      geom_bar() + #Setting the layer to a bar graph
      theme_classic() + #Setting the theme to a classic theme (no grid)
      scale_y_discrete(limits=c("Random 1", "Random 2", "Whole Brain", "Frontoparietal", "Salience Ventral Attention", "Default-Mode", "Dorsal Attention", "Limbic", "Somatomotor","Visual")) + #Setting up the columns to a specific order comes out. The order comes from most to least significant in the "Rest/Task" Pair. The order will be the same for all pairs.
      scale_fill_manual(values = c(fill_fail, fill_pass), guide = FALSE) + #Colors to use to fill the bar graph.
      labs(fill = "Fingerprinting result") + #Legend title
      theme(legend.position = "right") + #Legend position
      theme(axis.title.x=element_blank(), #Erasing the axis x title #Erasing the axis x ticks
            axis.title.y=element_blank()) #Erasing the title of the y axis.
    
  } else if (type == "between" & atlas == "Schaefer") {
    df1 <- df %>%
      select(id_camcan, contains(glue('non_fp_{modality1}_{modality2}_{type}'))) %>% #Selecting only the columns with "modality1/modality2"
      rename(
        "Frontoparietal" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_frontoparietal_{corr}"),
        "Salience Ventral Attention" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_salience_ventral_attention_{corr}"),
        "Default-Mode" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_default_mode_{corr}"),
        "Dorsal Attention" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_dorsal_attention_{corr}"),
        "Limbic" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_limbic_{corr}"),
        "Somatomotor" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_somatomotor_{corr}"),
        "Visual" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_visual_{corr}")
      ) %>% #Rename all columns to a easy acronym to avoid messing with levels instead
      pivot_longer(!id_camcan, names_to="network", values_to="identification") %>% #Pivoting the dataframe to a long format (ggplot won't work otherwise)
      mutate(label_pos = factor(case_when(identification == 1 ~ "Identified", identification == 0 ~ "Not identified", TRUE ~ NA_character_))) %>% #Adding a "factor" column for "identified/not_identified" for the "fill" of the the graph part ()
      na.omit() #Omitting rows with missing data
    
    #Define the graph to produce
    p <- df1 %>%
      ggplot(aes(y=network, fill=factor(label_pos, levels=c("Not identified", "Identified")))) + #Bar graph where 1 column = 1 network
      geom_bar() + #Setting the layer to a bar graph
      theme_classic() + #Setting the theme to a classic theme (no grid)
      scale_y_discrete(limits=c("Frontoparietal", "Salience Ventral Attention", "Default-Mode", "Dorsal Attention", "Limbic", "Somatomotor","Visual")) + #Setting up the columns to a specific order comes out. The order comes from most to least significant in the "Rest/Task" Pair. The order will be the same for all pairs.
      scale_fill_manual(values = c(fill_fail, fill_pass), guide = FALSE) + #Colors to use to fill the bar graph.
      labs(fill = "Fingerprinting result") + #Legend title
      theme(legend.position = "right") + #Legend position
      theme(axis.title.x=element_blank(), #Erasing the axis x title #Erasing the axis x ticks
            axis.title.y=element_blank()) #Erasing the title of the y axis.
    
  } else if (type == "within" & atlas == "Power") {
      df1 <- df %>%
          select(id_camcan, contains(glue('non_fp_{modality1}_{modality2}_{type}'))) %>%
          rename(
              "Uncertain" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_uncertain_{corr}"),
              "Salience" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_salience_{corr}"),
              "Cingular/Opercular" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_cingulo_opercular_{corr}"),
              "Dorsal attention" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_dorsal_attention_{corr}"),
              "Frontoparietal" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_frontoparietal_{corr}"),
              "Auditory" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_auditory_{corr}"),
              "Default-mode" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_default_mode_{corr}"),
              "Visual" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_visual_{corr}"),
              "Cerebellar" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_cerebellar_{corr}"),
              "Ventral attention" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_ventral_attention_{corr}"),
              "Sensory/Somatomotor" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_sensory_somatomotor_{corr}"),
              "Whole brain" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_Whole_brain_{corr}"),
              "Subcortical" = glue("non_fp_{modality1}_{modality2}_{type}_{atlas}_subcortical_{corr}")
          ) %>% #Rename all columns to a easy acronym to avoid messing with levels instead
          pivot_longer(!id_camcan, names_to="network", values_to="identification") %>% #Pivoting the dataframe to a long format (ggplot won't work otherwise)
          mutate(label_pos = factor(case_when(identification == 1 ~ "Identified", identification == 0 ~ "Not identified", TRUE ~ NA_character_))) %>% #Adding a "factor" column for "identified/not_identified" for the "fill" of the the graph part ()
          na.omit() #Omitting rows with missing data
      
      #Define the graph to produce
      p <- df1 %>%
          ggplot(aes(y=network, fill=factor(label_pos, levels=c("Not identified", "Identified")))) + #Bar graph where 1 column = 1 network
          geom_bar() + #Setting the layer to a bar graph
          theme_classic() + #Setting the theme to a classic theme (no grid)
          scale_y_discrete(limits=c("Whole brain", "Frontoparietal", "Default-mode", "Salience", "Dorsal attention", "Ventral attention", "Cingular/Opercular", "Auditory", "Sensory/Somatomotor", "Visual", "Cerebellar", "Subcortical", "Uncertain")) + #Setting up the columns to a specific order comes out. The order comes from most to least significant in the "Rest/Task" Pair. The order will be the same for all pairs.
          scale_fill_manual(values = c(fill_fail, fill_pass), guide = FALSE) + #Colors to use to fill the bar graph.
          labs(fill = "Fingerprinting result") + #Legend title
          theme(legend.position = "right") + #Legend position
          theme(axis.title.x=element_blank(), #Erasing the axis x title #Erasing the axis x ticks
                axis.title.y=element_blank()) #Erasing the title of the y axis.
      
  }
  
  return(p)
  
}
#---------------

# Part 1 - Calculating the percentage of people passing the fingerprint

### Calculating the percentage for the pair
perc_calc_1 <- percentage_fingerprinted_calculator(fpc_data_clean, glue("{modality1}"), glue("{modality2}"), glue("{atlas}"), glue("{corr}"), glue("{type}"))

### Export the table to a file
write.csv(perc_calc_1, file = glue("{output}/Table 1 - Percentage of recognition for {modality1}_{modality2}_{atlas}_{corr}_{type}.csv"))


#############################
# Part 2 - Stacked bar chart of fingerprinting status by network

### Generate the plot
bar_plot_1 <- bar_graph_generator(fpc_data_clean, glue("{modality1}"), glue("{modality2}"), glue("{atlas}"), glue("{corr}"), glue("{type}"))

### Export the figure to a file
ggsave(glue("{output}/Figure 1 - Bar plot of identification rate {modality1} {modality2} {atlas} {corr} {type}.png"), 
       plot = bar_plot_1, width = 7, height = 4, units="in", dpi=320)















