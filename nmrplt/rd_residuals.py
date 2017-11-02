import matplotlib.pyplot as plt
import numpy as np

def rd_residuals(df, res=None, field=None, ax=None):
    """
    Plots the residuals of the fit (y(data) - y(fit)).
    
    INPUTS:
    df:  multi-indexed dataframe
    res:  residue number
    field:  0 or 1, corresponding to multi-index level 1
    ax:  subplotting axis to be used
    """
    
    fields = np.unique(df.index.get_level_values(1))
    colors = ['blue','red']

    if ax == None:
        #fig, ax = plt.subplots(1, 1, figsize = (10, 7.5));
        ax = subplot2grid((1, 1), (0, 0))
    if res == None:
        res = 152
    if field == None:
        field = 0

    x = df.loc[res, fields[field]].index
    y_dat = df.loc[res, fields[field]]['R2eff_fit']
    y_fit = df.loc[res, fields[field]]['R2eff_calc']
    resid = y_dat - y_fit
    ax.plot(x, np.zeros(len(resid)), '--', c='black', lw=1, zorder=0)
    ax.scatter(x, resid, c=colors[field], marker='D', s=50, 
                zorder=1, alpha=.85, label=str(fields[field]))
    return
