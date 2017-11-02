import matplotlib.pyplot as plt
import numpy as np

def rd_curve(df, res=None, field=None, ax=None, typ=None):
    """
    This function plots the RD data and fit points.
    
    INPUTS:
    df:  multi-indexed pandas dataframe
    res:  residue number
    field: 0 or 1
    ax:  subplotting axis to be used
    typ:  'data' or 'fit'
    """
    
    fields = np.unique(df.index.get_level_values(1))
    colors = ['#0009FF', '#FF0000', '#6E72FF', '#FF6E6E']
    alpha=1.
    
    if ax == None:
        fig, ax = plt.subplots(1, 1, figsize = (10, 7.5));
    if res == None:
        res = 152
    if field == None:
        field = 0
    if typ == None or typ == 'data':
        x = df.loc[res, fields[field]].index
        y = df.loc[res, fields[field]]['R2eff_fit'].values
        yerr = df.loc[res, fields[field]]['R2eff_error'].values
        col = colors[field]
        ax.errorbar(x, y, yerr=yerr, label=str(fields[field])+' Data', zorder=1, ms=10, 
                   mfc=col, ecolor=col, alpha=alpha, fmt='o', mew=2, capsize=0);    
    elif typ == 'Fit' or typ == 'fit':
        x, y = df.loc[res, fields[field]].index, df.loc[res, fields[field]]['R2eff_calc']
        ax.plot(x, y, label=str(fields[field])+' Fit', zorder=0, lw=5, c=colors[field+2]);
    return