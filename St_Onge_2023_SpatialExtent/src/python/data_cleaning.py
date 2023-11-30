""" Spatial Extent Project - Data Cleaning script

The goal of this script is to clean the ADNI data as downloaded from the LONI servers
in May 2022. We specifically downloaded and cleaned:
- UCBERKELEYAV1451_04_26_22.csv (TAU AV1451 data)
- UCBERKELEYAV45_04_26_22.csv (Amyloid AV45 data)
- UCBERKELEYFBB_04_26_22.csv (Amyloid FBB data)
- PTDEMOG.csv (Participants' demographic information)
- APOERES.csv (Participants' ApoE4 status)
- UWNPSYCHSUM_01_23_23_07Sep2023.csv (Participants' cognition, memory and executive function)
- ADSP_PHC_COGN_05_06_22.csv (Participants' cognition, visuospatial and language)


"""

# Script set-up

import pandas as pd
import numpy as np
from sihnpy import spatial_extent as spex

sourcedata = '~/Desktop/vlpp_projects/St_Onge_2023_SpatialExtent/data/sourcedata'
derivatives = '~/Desktop/vlpp_projects/St_Onge_2023_SpatialExtent/data/derivatives'

#Dictionary of composite regions we use throughout the script
comp_regions = {"braak1_suvr":["CTX_LH_ENTORHINAL_SUVR", "CTX_RH_ENTORHINAL_SUVR"],

    "braak3_suvr":["CTX_LH_AMYGDALA_SUVR", "CTX_RH_AMYGDALA_SUVR",
    "CTX_LH_FUSIFORM_SUVR", "CTX_RH_FUSIFORM_SUVR", 
    "CTX_LH_PARAHIPPOCAMPAL_SUVR", "CTX_RH_PARAHIPPOCAMPAL_SUVR",
    "CTX_LH_LINGUAL_SUVR", "CTX_RH_LINGUAL_SUVR"],

    "braak4_suvr":["CTX_LH_INFERIORTEMPORAL_SUVR", "CTX_RH_INFERIORTEMPORAL_SUVR",
    "CTX_LH_MIDDLETEMPORAL_SUVR", "CTX_RH_MIDDLETEMPORAL_SUVR",
    "CTX_LH_ISTHMUSCINGULATE_SUVR", "CTX_RH_ISTHMUSCINGULATE_SUVR",
    "CTX_LH_CAUDALANTERIORCINGULARE_SUVR", "CTX_RH_CAUDALANTERIORCINGULATE_SUVR",
    "CTX_LH_INSULA_SUVR", "CTX_RH_INSULA_SUVR",
    "CTX_LH_POSTERIORCINGULATE_SUVR", "CTX_RH_POSTERIORCINGULATE_SUVR",
    "CTX_LH_ROSTRALANTERIORCINGULATE_SUVR", "CTX_RH_ROSTRALANTERIORCINGULATE_SUVR"],

    "braak5_suvr":["CTX_LH_LATERALOCCIPITAL_SUVR", "CTX_RH_LATERALOCCIPITTAL_SUVR",
    "CTX_LH_INFERIORPARIETAL_SUVR", "CTX_RH_INFERIORPARIETAL_SUVR",
    "CTX_LH_SUPERIORTEMPORAL_SUVR", "CTX_RH_SUPERIORTEMPORAL_SUVR",
    "CTX_LH_BANKSSTS_SUVR", "CTX_RH_BANKSSTS_SUVR",
    "CTX_LH_PRECUNEUS_SUVR", "CTX_RH_PRECUNEUS_SUVR",
    "CTX_LH_PARSOPERCULARIS_SUVR", "CTX_RH_PARSOPERCULARIS_SUVR",
    "CTX_LH_PARSORBITALIS_SUVR", "CTX_RH_PARSORBITALIS_SUVR",
    "CTX_LH_PARSTRIANGULARIS_SUVR", "CTX_RH_PARSTRIANGULARIS_SUVR",
    "CTX_LH_FRONTALPOLE_SUVR", "CTX_RH_FRONTALPOLE_SUVR",
    "CTX_LH_CAUDALMIDDLEFRONTAL_SUVR", "CTX_RH_CAUDALMIDDLEFRONTAL_SUVR",
    "CTX_LH_LATERALORBITOFRONTAL_SUVR", "CTX_RH_LATERALORBITOFRONTAL_SUVR",
    "CTX_LH_MEDIALORBITOFRONTAL_SUVR", "CTX_RH_MEDIALORBITOFRONTAL_SUVR",
    "CTX_LH_ROSTRALMIDDLEFRONTAL_SUVR", "CTX_RH_ROSTRALMIDDLEFRONTAL_SUVR",
    "CTX_LH_SUPERIORFRONTAL_SUVR", "CTX_RH_SUPERIORFRONTAL_SUVR",
    "CTX_LH_SUPERIORPARIETAL_SUVR", "CTX_RH_SUPERIORPARIETAL_SUVR",
    "CTX_LH_SUPRAMARGINAL_SUVR", "CTX_RH_SUPRAMARGINAL_SUVR",
    "CTX_LH_TRANSVERSETEMPORAL_SUVR", "CTX_RH_TRANSVERSETEMPORAL_SUVR"],

    "braak6_suvr":["CTX_LH_PERICALCARINE_SUVR", "CTX_RH_PERICALCARINE_SUVR",
    "CTX_LH_CUNEUS_SUVR", "CTX_RH_CUNEUS_SUVR",
    "CTX_LH_PARACENTRAL_SUVR", "CTX_RH_PARACENTRAL_SUVR",
    "CTX_LH_POSTCENTRAL_SUVR", "CTX_RH_POSTCENTRAL_SUVR",
    "CTX_LH_PRECENTRAL_SUVR", "CTX_RH_PRECENTRAL_SUVR"],
    
    "meta_roi_tau_jack":["CTX_LH_AMYGDALA_SUVR", "CTX_RH_AMYGDALA_SUVR",
    "CTX_LH_ENTORHINAL_SUVR", "CTX_RH_ENTORHINAL_SUVR",
    "CTX_LH_PARAHIPPOCAMPAL_SUVR","CTX_RH_PARAHIPPOCAMPAL_SUVR"
    "CTX_LH_FUSIFORM_SUVR", "CTX_RH_FUSIFORM_SUVR",
    "CTX_LH_INFERIORTEMPORAL_SUVR","CTX_RH_INFERIORTEMPORAL_SUVR",
    "CTX_LH_MIDDLETEMPORAL_SUVR","CTX_RH_MIDDLETEMPORAL_SUVR"],

    "meta_roi_ab_ozlen":["CTX_LH_ROSTRALANTERIORCINGULATE_SUVR", "CTX_RH_ROSTRALANTERIORCINGULATE_SUVR",
    "CTX_LH_PRECUNEUS_SUVR", "CTX_RH_PRECUNEUS_SUVR",
    "CTX_LH_MEDIALORBITOFRONTAL_SUVR", "CTX_RH_MEDIALORBITOFRONTAL_SUVR",
    "CTX_LH_ROSTRALMIDDLEFRONTAL_SUVR", "CTX_RH_ROSTRALMIDDLEFRONTAL_SUVR",
    "CTX_LH_INFERIORPARIETAL_SUVR", "CTX_RH_INFERIORPARIETAL_SUVR",
    "CTX_LH_SUPERIORFRONTAL_SUVR", "CTX_RH_SUPERIORFRONTAL_SUVR",
    "CTX_LH_POSTERIORCINGULATE_SUVR", "CTX_RH_POSTERIORCINGULATE_SUVR"]
    }

####################################################
# Data cleaning individual data sources

########################
## PET data
""" There are 3 spreadsheets to import: AV1451 (Tau), AV45 (Amyloid) and FBB (Amyloid).

To clean the data, we need to:
1) Import all three spreadsheets
2) Do some basic clean-up and variable prep
3) Normalize the data by the region of interest
"""

dict_pet_data = {} # To store the PET dataframes
list_pet_files = ['UCBERKELEYAV1451_04_26_22.csv', #List of files to import
                'UCBERKELEYAV45_04_26_22.csv', 
                'UCBERKELEYFBB_04_26_22.csv']

### 1. - Import PET data
#Since we need to repeat most of the same operations across, we can loop and store in a dictionary
for file in list_pet_files:
    try:
        data_pet = pd.read_csv(f'{sourcedata}/{file}') #Import file

        #Quick check that the data in the spreadsheet is not duplicated for each examdate
        if data_pet.duplicated(subset=['RID', 'EXAMDATE']).sum() > 0:
            print(f'ERROR: Duplicated visits for some participants in file {file}')

        #Store the data in the dictionary, with a specific key for each file
        if file == 'UCBERKELEYAV1451_04_26_22.csv':
            dict_pet_data['tau'] = data_pet
            print(f'AV-1451 data available for {data_pet["RID"].nunique()} participants')
        elif file == 'UCBERKELEYAV45_04_26_22.csv':
            dict_pet_data['av45'] = data_pet
            print(f'AV-45 data available for {data_pet["RID"].nunique()} participants')
        else:
            dict_pet_data['fbb'] = data_pet
            print(f'FBB data available for {data_pet["RID"].nunique()} participants')

    except OSError:
        print(f'No file to read. Skipping {file}')

### 2. - Initial clean-up of PET data

dict_cort_data = {} #To store the clean data in each brain region
dict_ref_data = {} #To store the reference PET info (positivity, ref region, etc.)
dict_comp_data = {} #To store the composite data (average calculated in specific regions)

# Let's clean each of the PET data
for tracer, pet_data in dict_pet_data.items():
    pet_data_cp = pet_data.copy() #Create a copy to modify
    
    #First, create a new column to count how many visits each participant has.
    # We do this by grouping participants, creating a cumulative count and finally adding 1
    # because Python starts at 0.
    pet_data_cp[f'{tracer}_pet_number'] = pet_data_cp.groupby('RID').cumcount().add(1)

    #Rename columns and set index to the RID and PET number label
    pet_data_clean = (pet_data_cp
        .rename(columns={'EXAMDATE':f'EXAMDATE_{tracer}',
                    #Amygdala is not cortical, but code is simplified by a lot if we write that
                        'LEFT_AMYGDALA_SUVR':'CTX_LH_AMYGDALA_SUVR', 
                        'RIGHT_AMYGDALA_SUVR':'CTX_RH_AMYGDALA_SUVR'})
        .set_index(['RID', f'EXAMDATE_{tracer}', f'{tracer}_pet_number'])
    )

    #Next, we need to isolate the reference regions for each tracer
    if tracer == 'tau':
        ref_region = pet_data_clean['INFERIORCEREBELLUM_SUVR']
    if tracer in ['av45', 'fbb']:
        ref_region = pet_data_clean['WHOLECEREBELLUM_SUVR']

    #We isolate the SUVR cortical values and divide them by the SUVR in the reference region
    dict_cort_data[tracer] = (pet_data_clean
        #Grab all the regions (should also grab the amygdala since we changed it before)
        .filter(regex='^CTX_.+SUVR', axis=1)
        #Drop the "unknown" freesurfer regions
        .drop(labels=['CTX_LH_UNKNOWN_SUVR', 'CTX_RH_UNKNOWN_SUVR'], axis=1)
        #Divide all values by the ref region
        .div(ref_region, axis=0)
    )

    #We then compute the reference regions we want for each tracer.
    if tracer == 'tau':
        dict_ref_data[tracer] = (pet_data_clean
            .filter(items=['INFERIORCEREBELLUM_SUVR', ''])
            .assign(tracer = 'FTP'))
    elif tracer == 'av45':
        dict_ref_data[tracer] = (pet_data_clean
            #Keep only the important variables
            .filter(items=['WHOLECEREBELLUM_SUVR', #Ref region SUVR
                        'SUMMARYSUVR_WHOLECEREBNORM', #Global SUVR
                        'SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF'], axis=1) #Positivity AV45
            #Rename the positivity to something clearer
            .rename(columns={"SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF":"AB_POSITIVITY_WCR"})
            .assign(tracer='AV45') #Assign a variable to keep track of tracer
            #Compute the centiloid value for this tracer
            .assign(centiloid = lambda x: (188.22 * x['SUMMARYSUVR_WHOLECEREBNORM']) - 189.16))
    elif tracer == 'fbb':
        dict_ref_data[tracer] = (pet_data_clean
            #Keep only important variables
            .filter(items=['WHOLECEREBELLUM_SUVR', #Ref region SUVR
                        'SUMMARYSUVR_WHOLECEREBNORM', #Global SUVR
                        'SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF'], axis=1) #Positivity FBB
            #Rename the positivity to something clearer and matching AV45
            .rename(columns={"SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF":"AB_POSITIVITY_WCR"})
            .assign(tracer='AV45') #Assign a variable to keep track of tracer
            #Compute the centiloid value for this tracer
            .assign(centiloid = lambda x: (157.15 * x['SUMMARYSUVR_WHOLECEREBNORM']) - 151.87))

    #Finally, we compute the composite data for the participants
    dict_comp_data[tracer] = pd.DataFrame(index=pet_data_clean.index)

    #For each composite region...
    for comp_name, regions in comp_regions.items():
        dict_comp_data[tracer][comp_name] = (pet_data_clean
            #Keep only the regions in the specific composite and average
            .filter(items=regions, axis=1)
            .mean(axis=1))

# Initial clean-up of the PET data is now done. We will get back to it later when the rest is also cleaned.

########################
## Clinical data

### 1. - Import clinical data

data_clin = pd.read_csv(f"{sourcedata}/DXSUM_PDXCONV_ADNIALL.csv")

### 2. - Clean clinical data

data_clin_clean = (data_clin
        .dropna(subset=['DIAGNOSIS']) #Drop missing
        .drop_duplicates(subset=['RID', 'EXAMDATE']) #Drop duplicate visits
        .filter(items=['RID', 'DIAGNOSIS', 'EXAMDATE'], axis=1) #Keep only the main variables
        .rename(columns={'EXAMDATE':'EXAMDATE_CLIN'}) #Rename the Examdate for clarity
        .assign(clin_number = lambda x: x.groupby('RID').cumcount().add(1))) 

########################
## Demographics data

### 1. - Import demographics data

data_demog = pd.read_csv(f"{sourcedata}/PTDEMOG.csv")

### 2. - Clean demographics data

data_demog_clean = (data_demog
    .filter(items=['RID', 'PTGENDER', 'PTDOBMM', 'PTDOBYY', 'PTEDUCAT'])
    .drop_duplicates(subset=['RID'])
    .assign(PTDOB = lambda x: 
            x['PTDOBYY'].astype(int).astype(str) + "-" +
            x['PTDOBMM'].astype(int).astype(str) + "-" +
            #The exact day of birth is not given to preserve identity.
            #We create a random day of birth to calculate the age later on
            f"{np.random.randint(1, 28)}")
    .drop(labels=['PTDOBMM', 'PTDOBYY'], axis=1))

########################
## ApoE4 data

### 1. - Import apoe data

data_apoe = pd.read_csv(f"{sourcedata}/APOERES.csv")

### 2. - Clean apoe data

data_apoe_clean = (data_apoe
    #In ADNI, ApoE status is in 2 variables
    .filter(items=['RID','APGEN1', 'APGEN2'], axis=1)
    #Create a binary variable if at least 1 apoe allele
    .assign(apoe_bin = lambda x: np.where((x['APGEN1'] == 4) | (x['APGEN2'] == 4), 'E4', "non_E4"))
    #Create a variable for full genotype
    .assign(apoe_genotype = lambda x: x['APGEN1'].astype(str) + "_" + x['APGEN2'].astype(str))
    #To avoid extra categories, reorder genotypes that are in another order
    .replace(to_replace={'apoe_genotype':{"4_3":"3_4", "3_2":"2_3", "4_2":"2_4"}})
    .drop(labels=['APGEN1', 'APGEN2'], axis=1)) #We don't need them anymore

########################
## Cognition data
#This one is slightly more tricky. There are two dataframes to merge:
## One with the memory and EF data, one with VS and LAN.

### 1. - Import cog data

data_cog_memef = pd.read_csv(f'{sourcedata}/UWNPSYCHSUM_12_13_21.csv')
data_cog_vislan = pd.read_csv(f'{sourcedata}/ADSP_PHC_COGN_05_06_22.csv')

### 2. Clean cognition data

data_cog_memef_clean = (data_cog_memef
        .filter(items=['RID', 'examdate', 'ADNI_MEM', 'ADNI_EF'], axis=1)
        .rename(columns={'examdate':'EXAMDATE_COG'}))

data_cog_vislan_clean = (data_cog_vislan
        .filter(items=['RID', 'EXAMDATE', 'PHC_LAN', 'PHC_VSP'], axis=1)
        .rename(columns={'EXAMDATE':'EXAMDATE_COG',
                        'PHC_LAN':'ADNI_LAN', 'PHC_VSP':'ADNI_VSP'}))

data_cog_clean = (data_cog_memef_clean.set_index(['RID', 'EXAMDATE_COG'])
        .merge(data_cog_vislan_clean.set_index(['RID', 'EXAMDATE_COG']), 
                    left_index=True, right_index=True, how='outer')
        .reset_index()
        .assign(cog_number= lambda x: x.groupby('RID').cumcount().add(1)))

########################
## Final inclusion 
# Here, we need to prepare the final data included in the study. We have few criteria,
# but they need to be set:
# 1) Must have at least 1 baseline tau timepoint
# 2) Must have amyloid timepoint within 2 years of the tau
# 3) Must have clinical data available within 2 years of the tau

### 1. - Merge amyloid in a single spreadsheet

data_ab_ref_merged = pd.concat( #We stack the two spreadsheets
    #A bit of a pain, but to rename the columns here, I need to copy the object first.
    [dict_ref_data["av45"].copy().reset_index().rename(columns={"EXAMDATE_av45":"EXAMDATE_AB",
                                            'av45_pet_number':'ab_pet_number'}),
    dict_ref_data['fbb'].copy().reset_index().rename(columns={'EXAMDATE_fbb':'EXAMDATE_AB',
                                            'fbb_pet_number':'ab_pet_number'})],
    axis=0, ignore_index=True
)

#Check if there are FBB/AV45 duplicates, and if all the data is not missing
#print(data_ab_ref_merged.reset_index().duplicated(subset=['RID','EXAMDATE_AB']).sum())
#print(data_ab_ref_merged['EXAMDATE_AB'].isnull().sum())

### 2. - Find the baseline tau data, and match the closest amyloid

tau_bls = (dict_ref_data['tau'].reset_index()
        .filter(items=['RID', 'EXAMDATE_tau', 'tau_pet_number']) #Keep only essential info
        .loc[lambda x: x['tau_pet_number'] == 1] #Restrict to baseline visits
        .drop(labels='tau_pet_number', axis=1) #No longer useful after restriction
        #Merge the amyloid data using a left merge, restrict to only date first
        .merge(data_ab_ref_merged.filter(items=['RID','EXAMDATE_AB'], axis=1),
            left_on='RID', right_on='RID', how='left')
        .dropna(subset=['EXAMDATE_AB'], axis=0) #If missing, can't compute a distance, so remove
        #Calculate distance in days between tau and amyloid scan
        .assign(tau_amyl_dist = lambda x: 
                ((pd.to_datetime(x['EXAMDATE_tau']) - pd.to_datetime(x['EXAMDATE_AB']))
                / np.timedelta64(1, 'D')).abs())
        #Isolate the closest amyloid scan to the baseline tau
        .loc[lambda x: x.groupby('RID')['tau_amyl_dist'].idxmin()]
        .loc[lambda x: x['tau_amyl_dist'] < (365*2)] #Remove participants where closest is >2 years
        #Next, we repeat this process with clinical diagnosis
        .merge(data_clin_clean, left_on='RID', right_on='RID', how='left') #Merge with data
        .dropna(subset=['EXAMDATE_CLIN'], axis=0) #Remove missing clinical dates
        #Calculate the distance in days between tau and clinical visit
        .assign(tau_clin_dist = lambda x:
                ((pd.to_datetime(x['EXAMDATE_tau']) - pd.to_datetime(x['EXAMDATE_CLIN']))
                / np.timedelta64(1, 'D')).abs())
        .loc[lambda x: x.groupby('RID')['tau_clin_dist'].idxmin()] #Find closest diagnosis
        .loc[lambda x: x['tau_clin_dist'] < (365*2)] #Remove visits more than 2 years away
)

### 3. Compute the spatial extent
""" Now that we have the final sample size, we can execute the steps from the spatial extent

"""

#### 3.1 - Preliminary cleaning
# Need to isolate the SUVR for the final IDS included

suvr_spex = (tau_bls.set_index('RID')
        .filter(items=[], axis=1) #Remove all columns except IDs
        #Isolate the baseline tau, and remove the non-cortical values, before merging.
        .merge(dict_cort_data['tau'].reset_index().loc[lambda x: x['tau_pet_number'] == 1]
            .drop(labels=['EXAMDATE_tau', 'tau_pet_number'], axis=1),
            left_on='RID', right_on='RID', how='left')
        #Set index so sihnpy ignores it
        .set_index('RID'))

#### 3.2 - sihnpy spatial extent - Threshold derivation step

#Estimate the GMM (2 components)
gm_estimations, clean_data = spex.gmm_estimation(data_to_estimate=suvr_spex, fix=False)
#Compute measures of each cluster and clean data
final_data, final_gm_dict, gmm_measures = spex.gmm_measures(cleaned_data=clean_data,
                                                            gm_objects=gm_estimations, fix=False)
#Computes abnormality probability for each sample
probability_data = spex.gmm_probs(final_data=final_data, 
                                final_gm_estimations=final_gm_dict, fix=False)
# (Optional) Outputs histograms of the densities for each cluster
#dict_figures = spex.gmm_histograms(final_data=final_data, 
#                               gmm_measures=gmm_measures, 
#                               probs_df=probability_data, type="density") 
#Derive thresholds based on probability.
thresh_df = spex.gmm_threshold_deriv(final_data=final_data,
                                    probs_df=probability_data, prob_threshs=[0.5], improb=1.0) 
#Don't forget to set an "improbable" value based on your data, if applicable.

#### 3.3 - sihnpy spatial extent - Threshold application step

#Derive the longitudinal tau data on which to apply the data
suvr_spex_long = (tau_bls
        .filter(items=['RID'], axis=1) #Remove all columns except IDs
        #Isolate the baseline tau, and remove the non-cortical values, before merging.
        .merge(dict_cort_data['tau'].reset_index(),
            left_on=['RID'], right_on=['RID'], how='left')
        #Set index so sihnpy ignores it
        .set_index(['RID', 'EXAMDATE_tau', 'tau_pet_number']))

#Apply the thresholds
data_to_apply_clean, thresh_data_clean = spex.apply_clean(data_to_apply=suvr_spex_long,
                        thresh_data=thresh_df) #Basic clean-up and reordering
dict_masks = spex.apply_masks(data_to_apply_clean=data_to_apply_clean,
                        thresh_data_clean=thresh_data_clean) #Create binary masks from thresholds
spex_metrics = spex.apply_index(data_to_apply_clean=data_to_apply_clean, 
                        dict_masks=dict_masks) #Computes spatial extent index

### 4. - Derive a spatial extent from a different threshold measure (2SD from CU AB-)

#### 4.1 - Need to identify the people in question

thresh_df_cu_neg = (tau_bls
        #Keep the amyloid visit linked to the initial tau visit and the diagnosis
        .filter(items=['RID', 'EXAMDATE_AB', 'DIAGNOSIS'], axis=1)
        #We're missing the positivity info to filter. Adding here.
        .merge(data_ab_ref_merged[['RID', 'EXAMDATE_AB', 'AB_POSITIVITY_WCR']], 
            left_on=['RID', 'EXAMDATE_AB'], right_on=['RID', 'EXAMDATE_AB'], how='left')
        .loc[lambda x: x['DIAGNOSIS'] == 1] #Keeping only CU
        .loc[lambda x: x['AB_POSITIVITY_WCR'] == 0] #Keeping only AB-
        .filter(items=['RID'], axis=1) #Keep only the IDs, we don't need the rest after filtering
        #Merge the SUVR data, keeping only the baseline tau
        .merge(dict_cort_data['tau'].reset_index().loc[lambda x: x['tau_pet_number'] == 1],
            left_on=['RID'], right_on=['RID'], how='left')
        #Remove the extra variables we don't need
        .drop(labels=['EXAMDATE_tau', 'tau_pet_number'], axis=1)
        .set_index('RID') #Set index to RID
        #Now, we need to compute the values we need, which is mean + 2SD in each column
        #We can abbreviate that with a lambda function.
        .agg([lambda x: x.mean() + ((x.std()) * 2)], axis=0)
        .melt() #Flips the data to long format, where each row is a brain region
        .rename(columns={'variable':'region', 'value':'thresh_2sd_cu_ab_neg'})
        .set_index('region') #Set the index to regions
)

#### 4.2 - Use sihnpy's spex apply functions to whole data
data_to_apply_clean_cu_neg, thresh_data_clean_cu_neg = spex.apply_clean(
                        data_to_apply=suvr_spex_long,
                        thresh_data=thresh_df_cu_neg) #Basic clean-up and reordering
dict_masks_cu_neg = spex.apply_masks(
                        data_to_apply_clean=data_to_apply_clean_cu_neg,
                        thresh_data_clean=thresh_data_clean_cu_neg) 
                        #Create binary masks from thresholds
spex_metrics_cu_neg = spex.apply_index(data_to_apply_clean=data_to_apply_clean_cu_neg, 
                        dict_masks=dict_masks_cu_neg) #Computes spatial extent index


### 5. - Merge the rest of the demographic data to the baseline tau data we just created

baseline_data = (tau_bls
        #First, we need to merge the cognition, finding the closest visit to the tau
        # No need to remove people over 2 years though
        .merge(data_cog_clean, left_on='RID', right_on='RID', how='left')
        .assign(tau_cog_dist = lambda x:
                ((pd.to_datetime(x['EXAMDATE_tau']) - pd.to_datetime(x['EXAMDATE_COG']))
                / np.timedelta64(1, 'D')).abs())
        .assign(tau_cog_dist = lambda x:
                np.where(x['tau_cog_dist'].notnull(), x['tau_cog_dist'], 9999))
        #Isolate closest cognition visit to the tau
        .loc[lambda x: x.groupby('RID')['tau_cog_dist'].idxmin()]
        #Turn back the 9999 values to missing for cognition
        .replace(to_replace={'tau_cog_dist':{9999:np.NaN}})
        #We can then merge the amyloid information, tau information,
        #  demographic and the apoe4 information
        #Let's start with demographic
        .merge(data_demog_clean, left_on='RID', right_on='RID', how='left')
        #We can now calculate the age, based on the age at the tau-PET
        .assign(PTAGE = lambda x:
                ((pd.to_datetime(x['EXAMDATE_tau']) - pd.to_datetime(x['PTDOB']))
                / np.timedelta64(1, 'Y')))
        .drop(labels=['PTDOB'], axis=1) #Don't need the birthday anymore
        #Next, we can merge the ApoE4 information
        .merge(data_apoe_clean, left_on='RID', right_on='RID', how='left')
        #Finally, we just need to merge back the amyloid reference (positivity + centiloid)
        #We need to merge on RID and EXAMDATE_AB (to avoid duplicate rows)
        .merge(data_ab_ref_merged, 
            left_on=['RID', 'EXAMDATE_AB'], right_on=['RID', 'EXAMDATE_AB'], how='left')
        #Merge the composite regions of tau
        .merge(dict_comp_data['tau'].add_suffix('_tau').reset_index()
                .loc[lambda x: x['tau_pet_number'] == 1]
                .drop(labels=['tau_pet_number'], axis=1), 
            left_on=['RID', 'EXAMDATE_tau'], right_on=['RID', 'EXAMDATE_tau'], how='left')
        #Merge the spatial extent at baseline
        .merge(spex_metrics.reset_index()
                .loc[lambda x: x['tau_pet_number'] == 1]
                .drop(labels=['tau_pet_number'], axis=1),
            left_on=['RID', 'EXAMDATE_tau'], right_on=['RID', 'EXAMDATE_tau'], how='left')
        #Merge the spatial extent at baseline from 2SD from CU AB-
        .merge(spex_metrics_cu_neg.reset_index()
                .loc[lambda x: x['tau_pet_number'] == 1]
                .drop(labels=['tau_pet_number'], axis=1),
            left_on=['RID', 'EXAMDATE_tau'], right_on=['RID', 'EXAMDATE_tau'], how='left')
        .set_index('RID')
)

### 6. - Prepare the longitudinal data
#Longitudinal data is slightly different to prepare. We just need to make sure that
# people in the longitudinal datasets are included in our list of participants

longitudinal_ab_data = (baseline_data.filter(items=[], axis=1)
    .merge(data_ab_ref_merged.set_index('RID'), left_index=True, right_index=True, how='left'))

longitudinal_taucomp_data = (baseline_data.filter(items=[], axis=1)
    .merge(dict_comp_data['tau'].add_suffix('_tau').reset_index().set_index('RID'),
        left_index=True, right_index=True, how='left'))

longitudinal_clin_data = (baseline_data.filter(items=[], axis=1)
    .merge(data_clin_clean.set_index('RID'), left_index=True, right_index=True, how='left'))

longitudinal_cog_data = (baseline_data.filter(items=[], axis=1)
    .merge(data_cog_clean.set_index('RID'), left_index=True, right_index=True, how='left'))

#print(longitudinal_ab_data.index.nunique()) #Check if unique RIDs is the same as planned
#print(longitudinal_taucomp_data.index.nunique()) #Check if unique RIDs is the same as planned
#print(longitudinal_clin_data.index.nunique()) #Check if unique RIDs is the same as planned
#print(longitudinal_cog_data.index.nunique()) #Check if unique RIDs is the same as planned

### 7. - Extra analyses 
""" Section containing extra data preparations run in addition to the original ones
"""

#### 7.1 - Apply spatial extent to the Braak stages

braak_suvr_spex = (tau_bls.set_index('RID')
        .filter(items=[], axis=1) #Remove all columns except IDs
        #Isolate the baseline tau, and remove the non-cortical values, before merging.
        .merge(dict_comp_data['tau'].reset_index().loc[lambda x: x['tau_pet_number'] == 1]
            .drop(labels=['EXAMDATE_tau', 'tau_pet_number', 'meta_roi_ab_ozlen'], axis=1),
            left_on='RID', right_on='RID', how='left')
        #Set index so sihnpy ignores it
        .set_index('RID'))

#### 7.2 - sihnpy spatial extent - Threshold derivation step

#Estimate the GMM (2 components)
gm_estimations_braak, clean_data_braak = spex.gmm_estimation(data_to_estimate=braak_suvr_spex,
                                                            fix=False)
#Compute measures of each cluster and clean data
final_data_braak, final_gm_dict_braak, gmm_measures_braak = spex\
    .gmm_measures(cleaned_data=clean_data_braak, gm_objects=gm_estimations_braak, fix=False)
#Computes abnormality probability for each sample
probability_data_braak = spex.gmm_probs(final_data=final_data_braak, 
                                final_gm_estimations=final_gm_dict_braak, fix=False)
# (Optional) Outputs histograms of the densities for each cluster
#dict_figures_braak = spex.gmm_histograms(final_data=final_data_braak, 
#                               gmm_measures=gmm_measures_braak, 
#                               probs_df=probability_data_braak, type="density") 
#Derive thresholds based on probability.
thresh_df_braak = spex.gmm_threshold_deriv(final_data=final_data_braak,
                                probs_df=probability_data_braak, prob_threshs=[0.5], improb=1.0) 
#Don't forget to set an "improbable" value based on your data, if applicable.

#### 7.3 - sihnpy spatial extent - Threshold application step

#Derive the longitudinal tau data on which to apply the data
braak_suvr_spex_long = (tau_bls
        .filter(items=['RID'], axis=1) #Remove all columns except IDs
        #Isolate the baseline tau, and remove the non-cortical values, before merging.
        .merge(dict_comp_data['tau'].reset_index(),
            left_on=['RID'], right_on=['RID'], how='left')
        #Set index so sihnpy ignores it
        .set_index(['RID', 'EXAMDATE_tau', 'tau_pet_number']))

#Apply the thresholds
data_to_apply_clean_braak, thresh_data_clean_braak = spex\
    .apply_clean(data_to_apply=braak_suvr_spex_long, thresh_data=thresh_df_braak) 
    #Basic clean-up and reordering
dict_masks_braak = spex.apply_masks(data_to_apply_clean=data_to_apply_clean_braak,
                        thresh_data_clean=thresh_data_clean_braak) 
                        #Create binary masks from thresholds
spex_metrics_braak = spex.apply_index(data_to_apply_clean=data_to_apply_clean_braak, 
                        dict_masks=dict_masks_braak) #Computes spatial extent index

#### 7.4 - Repeat the operations for the CU 2SD threshold

thresh_df_cu_neg_braak = (tau_bls
        #Keep the amyloid visit linked to the initial tau visit and the diagnosis
        .filter(items=['RID', 'EXAMDATE_AB', 'DIAGNOSIS'], axis=1)
        #We're missing the positivity info to filter. Adding here.
        .merge(data_ab_ref_merged[['RID', 'EXAMDATE_AB', 'AB_POSITIVITY_WCR']], 
            left_on=['RID', 'EXAMDATE_AB'], right_on=['RID', 'EXAMDATE_AB'], how='left')
        .loc[lambda x: x['DIAGNOSIS'] == 1] #Keeping only CU
        .loc[lambda x: x['AB_POSITIVITY_WCR'] == 0] #Keeping only AB-
        .filter(items=['RID'], axis=1) #Keep only the IDs, we don't need the rest after filtering
        #Merge the SUVR data, keeping only the baseline tau
        .merge(dict_comp_data['tau'].reset_index().loc[lambda x: x['tau_pet_number'] == 1],
            left_on=['RID'], right_on=['RID'], how='left')
        #Remove the extra variables we don't need
        .drop(labels=['EXAMDATE_tau', 'tau_pet_number'], axis=1)
        .set_index('RID') #Set index to RID
        #Now, we need to compute the values we need, which is mean + 2SD in each column
        #We can abbreviate that with a lambda function.
        .agg([lambda x: x.mean() + (x.std()) * 2], axis=0)
        .melt() #Flips the data to long format, where each row is a brain region
        .rename(columns={'variable':'region', 'value':'thresh_2sd_cu_ab_neg'})
        .set_index('region') #Set the index to regions
)

#### 4.4 - Repeat the threshold application for the Braak stages

#Apply the thresholds
data_to_apply_clean_braak_cu_neg, thresh_data_clean_braak_cu_neg = spex\
    .apply_clean(data_to_apply=braak_suvr_spex_long, thresh_data=thresh_df_cu_neg_braak) 
    #Basic clean-up and reordering
dict_masks_braak_cu_neg = spex.apply_masks(data_to_apply_clean=data_to_apply_clean_braak_cu_neg,
                        thresh_data_clean=thresh_data_clean_braak_cu_neg) 
                        #Create binary masks from thresholds
spex_metrics_braak_cu_neg = spex.apply_index(data_to_apply_clean=data_to_apply_clean_braak_cu_neg, 
                        dict_masks=dict_masks_braak_cu_neg) #Computes spatial extent index

### 8. - Output the results to file

baseline_data.to_csv(f'{derivatives}/baseline_data.csv')
spex_metrics.to_csv(f'{derivatives}/spex_metrics_gmm_long.csv')
spex_metrics_cu_neg.to_csv(f'{derivatives}/spex_metrics_cu_neg_long')
spex_metrics_braak.to_csv(f'{derivatives}/spex_metrics_gmm_braak_long.csv')
spex_metrics_braak_cu_neg.to_csv(f'{derivatives}/spex_metrics_cu_neg_long.csv')
suvr_spex_long.to_csv(f'{derivatives}/suvr_values_long.csv')
braak_suvr_spex_long.to_csv(f'{derivatives}/suvr_values_braak_long.csv')
longitudinal_ab_data.to_csv(f'{derivatives}/longitudinal_ab_data.csv')
longitudinal_taucomp_data.to_csv(f'{derivatives}/longitudinal_tau_comp_data.csv')
longitudinal_clin_data.to_csv(f'{derivatives}/longitudinal_clin_data.csv')
longitudinal_cog_data.to_csv(f'{derivatives}/longitudinal_cog_data.csv')
dict_masks['thresh_0.5'].to_csv(f'{derivatives}/long_spex_gmm_mask.csv')
dict_masks_cu_neg['thresh_2sd_cu_ab_neg'].to_csv(f'{derivatives}/long_spex_cu_neg_mask.csv')
dict_masks_braak['thresh_0.5'].to_csv(f'{derivatives}/long_spex_gmm_braak_mask.csv')
dict_masks_braak_cu_neg['thresh_2sd_cu_ab_neg']\
    .to_csv(f'{derivatives}/long_spex_cu_neg_braak_mask.csv')
thresh_df.to_csv(f'{derivatives}/regional_spex_thresholds.csv')
thresh_df_braak.to_csv(f'{derivatives}/braak_spex_thresholds.csv')
thresh_df_cu_neg.to_csv(f'{derivatives}/cu_neg_spex_thresholds.csv')
thresh_df_cu_neg_braak.to_csv(f"{derivatives}/braak_cuneg_spex_thresholds.csv")