import os
import re
import pandas as pd

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
