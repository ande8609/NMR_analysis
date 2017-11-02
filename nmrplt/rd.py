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


def rd_agg(df, res=None, ax1=None, ax2=None ):
    """
    Plots the RD curves and residuals together.
    
    INPUTS:
    df:  multi-indexed dataframe
    res:  residue number
    axs:  two subplotting axes to be used, one for curves
          one for residuals.
    """
    
    if res == None:
        res = 152

    if ax1 == None or ax2 == None:
        ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=3);
        ax2 = plt.subplot2grid((4, 1), (3, 0), sharex=ax1);
    
    xmax = 1300
    # First Axes: Data + Fit
    rd_curve(df, res, typ='data', ax=ax1, field=1)
    rd_curve(df, res, typ='fit', ax=ax1, field=1)
    rd_curve(df, res, typ='data', ax=ax1, field=0)
    rd_curve(df, res, typ='fit', ax=ax1, field=0)
    ax1.legend(framealpha = 1., loc=7);
    ax1.set_ylabel(r"$R_{2,eff} (\frac{1}{s})$", fontsize=20)
    ax1.set_title("Residue %s, RD Curves" %(res), fontsize=20)
    
    # Second Axes: Residuals
    rd_residuals(df, res, ax=ax2, field=0)
    rd_residuals(df, res, ax=ax2, field=1)
    ax2.set_xlabel(r"$\tau_{CPMG} (ms)$", fontsize=20);
    ax2.text(1100, 0.0, 'Residuals', transform=ax2.transData, fontsize=12)
    plt.xlim(-100, xmax);
    plt.tight_layout();
    return