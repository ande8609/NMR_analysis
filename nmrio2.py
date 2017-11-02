# PYTHON NMR I/O FUNCTIONS 
# Author:  Jeremy Anderson, 06/20/2016


import pandas as pd
import os
import fnmatch
import re
import matplotlib as mpl

###################
# General imports #
###################

def ccp_tsv_in(f):
    """
    Description:  Given a filepath, f, to a *.tsv peak list saved through ccpnmr imports
    the peak list and returns a pandas dataframe indexed by res number.
    
    Notes:  Add in a column of resnames.
    """
    
    df = pd.read_csv(f,sep='\t',usecols=[2,3,5,6,7,8,9],index_col=2)
    index = [int(re.findall(r'\d+',i)[0]) for i in df.index.values]
    res = [re.findall(r'[A-Za-z]{1,3}',i)[0] for i in df.index.values]
    df['AA'] = pd.Series(res)
    df.index = index
    return df.transpose().sort_index(axis=1)

def shift_import(path, isdir=True):
    """
    Description:  A function that imports a peak list or lists and returns the 1H and 15N 
    chemical shift lists.
    
    Input:  path - path to the peak list or directory
    isdir; boolean True/False, true if a directory is passed, false if not.
    """
    
    if isdir == False:
        big_df = ccp_tsv_in(path)
        h1_shift = big_df.xs('Position F1')
        n15_shift = big_df.xs('Position F2')
    else:
        big_df = read_tsv_peak_dir(path)
        h1_shift = big_df.xs('Position F1', level='feature')
        n15_shift = big_df.xs('Position F2', level='feature')
    return h1_shift, n15_shift

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

def hsqc_panel_in(img_dir):
    """
    Description:  Given a directory containing only panels of an HSQC spectrum, preferably
    six panels, this function imports each image using matplotlib and returns a list of 
    ndarrays.
    
    Input:  img_dir; image directory
    """
    
    
    # Image import: HSQC Panels
    hsqc_imgs = list()
    for f in os.listdir(img_dir):
        hsqc_imgs.append(mpl.image.imread(os.path.join(img_dir, f)))
        
    return hsqc_imgs

######################
# relax file imports #
######################

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

def relax_pardir(data_dir):
    """
    Description:  Given a relax analysis directory, imports all the parameters from the
    fit for every resn and returns a pandas dataframe indexed by resn with parameters and
    errors as col names.
    
    Inputs:
      data_dir:  relax analysis directory path.
    """
    
    import os, re
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
