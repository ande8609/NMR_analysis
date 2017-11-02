import os
import fnmatch
import re
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



def relax_pardir(data_dir):
    """
    Description:  Given a relax analysis directory, imports all the parameters from the
    fit for every resn and returns a pandas dataframe indexed by resn with parameters and
    errors as col names.
    
    Inputs:
      data_dir:  relax analysis directory path.
    """
    
    df_list = list()
    files = [f for f in os.listdir(data_dir) if re.search(r'[Aa-z0-9]{1,5}\.out', f)]
    for f in files:
        name = (os.path.splitext(os.path.basename(f))[0])
        temp_file = pd.read_table(os.path.join(data_dir, f), delim_whitespace=True, 
                                  skiprows=2, na_values='None', quotechar="'")
        header = list(temp_file.columns[1:])
        temp_head = [name+'_'+i for i in header[5:]]
        header[5:] = temp_head
        temp_file = temp_file.iloc[:,:-1]
        temp_file.columns = header
        df_list.append(temp_file)
    big_df = df_list[0]
    for n, f in enumerate(df_list):
        if n != 0:
            big_df = big_df.merge(f)
    big_df.set_index(header[1], inplace=True)
    return big_df.sort_index()
