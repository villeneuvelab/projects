import pandas as pd
import os
import numpy as np

#Once you have the windows data, the last part missing will just be to split the dataframe in windows (and store each dataframe in a dictionary with the right window).

#Import the fingerprinting file with the fp columns
fp_file = pd.read_csv("~/Desktop/Fingerprinting_Project/CAMCAN/05_Stats/01_Datasets/fpc_schaefer_pcor_data_clean.csv")

#We set the index and remove all columns but the identifiability
non_fp_only = fp_file.set_index("id_camcan")\
    .filter(like="non_fp")

#Flip to long format
long_fmt_fp = (non_fp_only
    .reset_index()
    .melt(id_vars=["id_camcan"]) #Flips and all columns become a single var
    .set_index('id_camcan') #Sets index
    #Change some variable names with extra "_" (otherwise column number is not equal)
    .replace(to_replace={"variable":{"Whole_brain":"wholebrain",
        "default_mode_limbic":"defaultmodelimbic",
        "default_mode":"defaultmode",
        "dorsal_attention":"dorsalattention",
        "salience_ventral_attention":"salienceventralattention"}},
        regex=True))

print(long_fmt_fp["variable"].value_counts())

#Split the 1 flipped variable in the different components
long_fmt_fp[["non", "fp", "mod1", "mod2", 
"type", "atlas", "network", "corr"]] = long_fmt_fp['variable']\
    .str.split("_", -1, expand=True)

#Remove the variables that are not informative
reduced_data = (long_fmt_fp
    .drop(labels=['variable', 'non', 'fp', 'atlas', 'corr'], axis=1)
    #Drop rows where the identifiability is missing. It's the only column where there can be stuff missing
    .dropna(how="any", axis=0)
)

########################### HERE WOULD BE THE WINDOWS OPERATION
# Import the data of all windows (like a for import loop, keep the ID column as a simple list, and store the list in a dictionary with key "window{i}")
dict_window_red = {}
for i, file in enumerate(os.listdir("~/Desktop/window_id_fpc/150_25")):
    #First, import each window, and store the first column as a list
    df = pd.read_csv(f"~/Desktop/window_id_fpc/150_25/{file}")
    list_window = df['id_camcan'].to_list()
    print(len(list_window))

    #Then, take subset the main dataframe so we only keep the data for each individual window.
    window_df = reduced_data.loc[list_window]
    print(len(window_df.index))

    #Finally, store the resulting window data in the dictionary
    dict_window_red[f'window{i + 1}'] = window_df
    #We force a plus 1 to facilitate the usage (instead of 0 indexed list)

###########################
#dict_window_red['window0'].reset_index()['id_camcan'].value_counts() #Seems to work. 100 unique IDs in the dataset

#Final step is to generate the percentage for each window.
# Basically a for loop on the dictionary we create and then we do the operation below
# And store each in a dictionary

dict_window_calc_perc = {}
for keys, items in dict_window_red.items():
    #For each subsetted data, we perform operations
    calc_perc = (items
        #Group by columns of interest (we want percentages for each unique group)
        .groupby(by=['mod1', 'mod2', 'type', 'network'])
        #Aggregate the data by counting the total cells and by summing the 1s
        .agg({"value":["count", "sum"]})
        #Agg creates a 2-level columns index. We drop the extra one
        .droplevel(axis=1, level=0)
        #Compute percentage for each group of interest, and store in a column.
        .assign(perc_rat = lambda x: (x['sum'] / x['count']))
        .assign(se_rat = lambda x: 1.96*np.sqrt((x['perc_rat'] * (1 - x['perc_rat'])) / x['count']))
        .assign(perc = lambda x: x['perc_rat'] * 100)
        .assign(se = lambda x: x['se_rat'] * 100)
        #Create final column to keep track
        .assign(window = keys))

    dict_window_calc_perc[keys] = calc_perc

dict_window_calc_perc['window1']
final_df = pd.concat(dict_window_calc_perc.values())

#Once the loop is done, we just need to concatenate everything together. We will then have all the data for all of the windows in long form, easy to use.
#Not sure if index will cause issue here. Worst case force reset in the loop

#We can then just export it for R

final_df.to_csv("~/Desktop/calc_perc_all_wind_150_25.csv")