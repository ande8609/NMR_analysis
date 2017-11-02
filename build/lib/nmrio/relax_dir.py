import os
import fnmatch
import pandas as pd

def relax_dir(directory):
    """
    Given a relax analysis directory containing 'disp*.out' files, this function reads all such
    files and returns a multi-indexed dataframe with the first level as residue number ('Residue') 
    the second level field strength ('Field') and the third index level cpmg pulsing frequency
    ('nu_cpmg').
    """
    
    res_list = list()
    df_list = list()
    for f in os.listdir(directory):
        if fnmatch.fnmatch(f, "disp*.out"):
            res_list.append(int(re.findall("[0-9]+",f)[0]))
            df_list.append( relax_file(directory, f) )
    
    big_df = pd.concat( df_list, axis=0, keys=res_list, names=['Residue', 'Field', 'nu_cpmg'] ).sort_index()
    # Needed for quick retrieval when slicing
    big_df.sortlevel(inplace=True)
    return big_df
