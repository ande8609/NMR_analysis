import os
import fnmatch
import pandas as pd
import re


def read_tsv_peak_dir(data_dir):
    """
    A function to read all peak list files of the form '*list.tsv' in the 
    given data directory.  This version also has functionality hard coded
    in it that pulls the spectrometer frequency out, which is predicated
    on the filename starting with that frequency.
    """
   
    # Finding all peak list files, *list.tsv, in the current directory
    # Also pulling out the 
    df_list = list()
    keys = list()
    for f in os.listdir(data_dir):
        if fnmatch.fnmatch(f, '*list.tsv'):
            df_list.append(ccp_tsv_in(os.path.join(data_dir, f)))
            keys.append(re.findall(r'\d+', f)[0])
    
    # Making a large dataframe from all peak list files
    big_df = pd.concat(df_list,axis=0,keys=keys,names=['peak list','feature'])
    return big_df
