"""Merging Fingerprinting output

The goal of the current package is to extract functional connectivity data from already scrubbed and cleaned fMRI data. The goal is simply to extract the data using specifical parcellations.

To do:
    - WARNING: Bug when running in variance data. For some reason, the index columns are not properly assigned and it generates duplicates in the final data. To fix.
        - Update: I'm an idiot. The reason why it duplicated is that the indices were not fixed in the variance script, to the 0 column is actually just numbers (not the ID)
        We fix this by allowing the user to input their own index, but the default is that it grabs the 0 column.
    - Customize to allow different parcellation options?
    - Create an option to select only certain subjects instead of a whole folder.
    - (Irrelevant) Create a timer for the extraction to know how much time it takes
    - Finish documenting the code.
    - (DONE) Make more general so you can merge the var dataframes
    - Probably could simplify the code by merging all the identifying string when building the class object. This would avoid always forgetting a "corr" or whatever else.
"""

import pandas as pd 
import numpy as np
import argparse
import os

def main():


    print('--------------------------------------------')
    print('-----Merging the dataframes with Python-----')
    print('--------------------------------------------')

    print('Parsing the arguments...')
    args = parsing()
    print('Verifying input...')
    verify_input(args)
    print('Generating an object of class fpc_df...')
    fpc_m = fpc_df(args.folder, args.output, args.type, args.index, args.corr, args.parcellation, args.extra, args.version)
    print('Importing the dataframes to merge...')
    fpc_m.import_df_names()
    print('Merging dataframes...')
    fpc_m.merge_dataframes()
    print('Exporting final dataframe...')
    fpc_m.export_dataframe()
    
    print('Done!')

    return None

def parsing():
    """
    The parsing() function needs the following parameters:
        sub_list            : Positional argument. Path leading to the folder with the subjects to process.
        output              : Positional argument. Path where the connectivity matrices should be output.
    """
    parser = argparse.ArgumentParser(description='Parsing the arguments for the matrix extraction')
    parser.add_argument('--folder', help='Path to dataframes to merge')
    parser.add_argument('--output', help='Path for output of the extraction')
    parser.add_argument('-t', '--type', default='within', choices=['within', 'between'], help='Whether to use "within" or "between" network connectivity for the FP.')
    parser.add_argument('-i', "--index", default=0, help="Index on which to merge all the datasets. By default, it will use the column in position 0 of each dataset. Technically, it should support any type of string that can be read in a python list (because that's what we feed it to), but this has to be tested.")
    parser.add_argument('-c', '--corr', help='Type of correlation ')
    parser.add_argument('-p', '--parcellation', help='Type of parcellation used')
    parser.add_argument('-e', '--extra', help='Extra string to add to the dataframe')
    parser.add_argument('-v', '--version', choices=['fpc', 'var'], help='Whether the script is used to merge fpc or var dataframes')

    args = parser.parse_args()

    return args

def verify_input(args):
    """
    Verify_input() does just that: it verifies that the arguments given to the parser make sense. The parser already validates the choice of connectivity measure and the choice of the parcellation.

    The current function simply verifies the paths given with output and sub_list (if they exist)
    """

    #Verify directories
    for paths in [args.folder, args.output]:
        if (os.path.isdir(paths) == True):
            pass
        else:
            raise SystemExit(f'The path given for the subject list ({paths}) does not exist. Please verify the path.')


    return None

class fpc_df:

    def __init__(self, folder, output, type, index, corr, parcellation, extra, version):
        self.folder = folder #Path where the files are
        self.output = output #Path for the output
        self.type = type #What type of measure we use for the fingerprinting
        self.index = index #What index column should we use?
        self.corr = corr #Type of correlation used in the fingerprinting. Effectively, just used as a string for the output name
        self.parce = parcellation #The parcellation used. Effectively, just used as a string for the output name.
        self.extra = extra #Any other string you need to add to the 
        self.ver = version #Whether the script should give a "fpc_" or "var_" name to the file
        
        return None
    
    def import_df_names(self):
        """
        This function scowers the folder given to the function, finds every .csv folder and import them as pandas dataframes. They are then added to the files list, which will be used in the next method. 
        """
        
        self.files = []

        for file in os.listdir(self.folder):
            file = os.path.join(self.folder, file)

            if (os.path.isfile(file) and file.endswith('.csv')):
                df = pd.read_csv(file, index_col=[self.index])
                self.files.append(df)
            else:
                print(f'Cannot process "{os.path.basename(file)}". Not a .csv file.')

        return None
    
    def merge_dataframes(self):
        """
        This function intakes the files instance of the fpc_df. It takes the first dataframe (list index 0) and sets it as the "primary dataframe". This is because the merge function of pandas needs a dataframe to merge on. We then recursively merge the dataframes in the files instance to the primary dataframe with an outer method (in case there are cases missing from one dataframe). All the dataframe have the same type of index, so we can just use the left_index/right_index arguments. 
        """

        self.len_files = len(self.files)
        self.prim_df = self.files[0]

        self.rest_df = self.files[1:self.len_files]

        for df in self.rest_df:
            self.prim_df = self.prim_df.merge(df, left_index=True, right_index=True, how='outer')

        return None

    def export_dataframe(self):
        """
        This method simply exports the final merged csv.
        """
        if self.extra:
            self.prim_df.to_csv(f'{self.output}/{self.ver}_{self.type}_{self.parce}_{self.corr}_{self.extra}.csv')
        else:
            self.prim_df.to_csv(f'{self.output}/{self.ver}_{self.type}_{self.parce}_{self.corr}.csv')

        return None


if __name__ == "__main__":
	main()