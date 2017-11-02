import pandas as pd
import re

def ccp_tsv_in(f, wonky=False):
    """
    Description:  Given a filepath, f, to a *.tsv peak list saved through ccpnmr imports
    the peak list and returns a pandas dataframe indexed by res number.
    
    Notes:  Add in a column of resnames.
    """
    if wonky == False:
        df = pd.read_csv(f,sep='\t',usecols=[2,3,5,6,7,8,9],index_col=2)
    else:
        df = pd.read_csv(f,sep='\t',usecols=[2,3,4,6,7,8,9],index_col=2)
    index = [int(re.findall(r'\d+',i)[0]) for i in df.index.values]
    res = [re.findall(r'[A-Za-z]{1,3}',i)[0] for i in df.index.values]
    df['AA'] = pd.Series(res)
    df.index = index
    return df.transpose().sort_index(axis=1)
