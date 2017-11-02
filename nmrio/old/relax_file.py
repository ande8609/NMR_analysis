import os
import pandas as pd

def relax_file(directory, filename):
    """
    A function that takes a directory and a 'disp*.out' file from a relax analysis and returns
    a dataframe multi-indexed by field strength at the first level and nu_cpmg at the second
    level.
    """
    df = pd.read_table( os.path.join(directory, filename), skiprows=1, usecols=range(2,7), header=None, 
                             names=['Field','nu_cpmg','R2eff_fit','R2eff_calc','R2eff_error'], 
                             delim_whitespace=True, na_values='-', index_col=[0,1] )
    return df
