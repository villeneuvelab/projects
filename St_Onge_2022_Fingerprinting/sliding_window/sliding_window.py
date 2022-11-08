""" Sliding-window module

Sliding-window analysis. The code is a python adaptation
from the code used by Vasa et al. (2018).

It outputs a spreadsheet for each window, where the IDs of
that window and the sorting variable are output. It also
outputs some descriptive stats for the sorting variable in
each window.
"""

import pandas as pd 
import math

def import_whole_sample(path_file, sorting_var):
    """Simple function importing the data with the participant ids and the variable to sort on.

    Parameters
    ----------
    path_file : str
        Path to the .csv file where the IDs and sorting variable are in.
    sorting_var : str
        Name (string) of the column to use for sorting.

    Returns
    -------
    pandas.DataFrame
        Dataframe containing the rows sorted by the variable of interest. Always sorted from
        smallest to biggest value.
    """

    df_whole_sample = pd.read_csv(path_file, index_col=0, usecols=[sorting_var])

    df_sorted = df_whole_sample.sort_values(by=[sorting_var], ascending=True, axis=0)

    return df_sorted

def compute_n_bins(df_sorted, ww, sts):
    """Function to compute how many bins need to be done, depending on the window size and step
    size.

    The window size of a window is the number of participant in a given window. The step size
    is how much we should move before selection the next window. The a large step size means a
    smaller overlap between windows, while a small step size meas a larger overlap between windows. 

    The original methodology and code was developped by Vasa et al. (2018). The only adaptation
    is that when window parameters do not divide a sample in a clean number of windows (i.e.,
    there is a remainder of participants), we add an extra window to add the rest of the
    participants.

    Parameters
    ----------
    df_sorted : pandas.DataFrame
        Dataframe of participants, sorted by the variable of interest.
    ww : int
        Integer representing the window size; how many participants should be in each window.
    sts : int
        Integer representing the step size; how many participants should be skipped in selecting
        participants in the next window.

    Returns
    -------
    int
        Returns the number of bins that the sample should be divided in.

    Raises
    ------
    SystemExit
        In the eventuality that the number of participants in the dataframe is less than
        the window size, the script will throw an error. It's a bit primitive, but is
        more meant to avoid typing errors.
    """
    num_part = len(df_sorted) #Total number of participants.

    if (num_part - ww) < 0:
        raise SystemExit('  ERROR: The number of participants is less than the window size.')

    print("Testing for clean division...")
    if ((num_part - ww) % sts > 0):
        print("    - Division is clean: no adjustment needed")
        nbin = math.ceil((num_part - ww)/sts)
        print(f'    - Number of windows: {nbin}')
    else:
        print("    - Doesn't divide cleanly. Last window may have a different number of " +
        "participants.")
        nbin = math.ceil((num_part - ww)/sts) + 1
        print(f"    - Number of windows: {nbin}")

    return nbin

def sliding_windows(df_sorted, ww, sts, nbin):
    """Actual function computing which participants should be in which windows. Based on the
    original function in R in Vasa et al. (2018).

    The main difference is that in R, indexing starts at 1, while it starts at 0 in Python.
    The code written accounts for that difference. We also had issues where the script would
    ignore the remainder of participants after the window. So we added a small clause where, if we
    reach the last window, all remaining participants are included.

    Parameters
    ----------
    df_sorted : pandas.DataFrame
        Dataframe with sorted IDs
    ww : int
        Integer representing the window size; how many participants should be in each window.
    sts : int
        Integer representing the step size; how many participants should be skipped in selecting
        participants in the next window.
    nbin : int
        Integer indicating how many windows we should be making.

    Returns
    -------
    dict
        Returns a dictionary, where the keys are the name the file we will export should bear,
        while the values are pandas.DataFrame with 1 column for ID and 1 column being the sorting
        variable.
    """

    w_store = {} #To store the windows
    for bin_num in range(0, nbin):
        bin_id = bin_num + 1
        print(f'    Creating window {bin_id}')

        if bin_id == nbin:
            #For some reason I couldn't figure out, the code would ignore the "remainder" 
            # participants in the last window. Here, I force the last window to take any remaining
            # participants so no one is left out of the analysis.
            bin_list = df_sorted.iloc[(sts*(bin_id-1)):]
        else:
            #In the original R code, the code below starts with 1+. However, Python is a 0 indexed
            # language, so it becomes unecessary here.
            bin_list = df_sorted.iloc[(sts*(bin_id-1)):(ww+sts*(bin_id-1))]

        #Here we just make a name for the file. We adapt a bit so the length of the names
        # are exactly equal in character length (easier for users to use regex later).
        if bin_id >= 10:
            bin_name = f"ww{ww}_sts{sts}_w{bin_id}"
        else:
            bin_name = f"ww{ww}_sts{sts}_w0{bin_id}"

        w_store[f'{bin_name}'] = bin_list

    return w_store

def window_metrics(w_store, sorting_var):
    """Summary function computing descriptive statistics for the variable used for sorting in each
    window.

    Each window is made by a different set of participants. It is often useful to get descriptive
    statistics on who is included in the different windows, either to report the composition of
    each window and/or to use on graphs as an indicative measure (instead of using "window no.1")

    Parameters
    ----------
    w_store : dict
        Dictionary of windows from the previous function, where the keys are the names of the files
        and the values are the actual IDs and sorting variable in a DataFrame.
    sorting_var : str
        String of the name of the column that was used to sort participants.

    Returns
    -------
    pandas.DataFrame
        Dataframe with mean, standard deviation, median, minimum and maximum values for the sorting
        variable for all windows with the current window parameters.
    """
    data_window_metrics = pd.DataFrame(columns=['window', 'measure', 'value'])

    for bin_name, window in w_store.items():
        tmp_df = pd.DataFrame(data={'window':[bin_name, bin_name, bin_name, bin_name, bin_name], 
            'measure':['mean', 'std', 'median', 'min', 'max'], 
            'value':[window[sorting_var].mean(), window[sorting_var].std(),
                window[sorting_var].median(), window[sorting_var].min(), 
                window[sorting_var].max()]})

        data_window_metrics = pd.concat([data_window_metrics, tmp_df])

    return data_window_metrics

def output_window(output, ww, sts, w_store, data_window_metrics):
    """Simple function exporting the IDs of each window in separate spreadsheets
    and exporting the computed measures within each window.

    Parameters
    ----------
    output : str
        String representing the path to the folder where the results should be output.
    w_store : dict
        Dictionary of file names and dataframes containing the IDs
    data_window_metrics : pandas.DataFrame
        Dataframe containing the descriptive statistics of each window
    """

    for bin_name, window in w_store.items():
        window.to_csv(f"{output}/{bin_name}.csv")

    data_window_metrics.to_csv(f"{output}/data_window_metrics_ww{ww}_sts{sts}.csv")
