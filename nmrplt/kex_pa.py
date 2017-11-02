import matplotlib.pyplot as plt
import numpy as np

def kex_pa(df, res=None, ax1=None, ax2=None, ax3=None, hist=True):
    """
    Description:  Plotting kex vs. pA CPMG parameters.
    
    Inputs:  
      df:  a pandas dataframe indexed by resn, columns named by ccpnmr.
      res:  residue number.
      ax1:  matplotlib axes object, where main data will be plotted.
      ax2:  matplotlib axes object, where kex histogram is plotted.
      ax3:  matplotlib axes object, where pA hist is plotted.
      hist:  boolean T/F, whether to plot hists or not.
    """
    
    colors = ['blue','red']

    if ax1 == None:
        #fig, ax = plt.subplots(1, 1, figsize = (10, 7.5));
        ax1 = plt.subplot2grid((4, 4), (0, 0), rowspan=3, colspan=3);
        ax2 = plt.subplot2grid((4, 4), (3, 0), colspan=3);
        ax3 = plt.subplot2grid((4, 4), (0, 3), rowspan=3);
        #plt.figure(figsize=(10, 7.5));
    if res == None:
        res = 152

    # Selecting data
    x = df['kex_value'].values
    y = df['pA_value'].values
    xerr = df['kex_error'].values
    yerr = df['pA_error'].values
    a = df.loc[res,'kex_value':'kex_error']
    b = df.loc[res,'pA_value':'pA_error']
    # Plotting all data
    ax1.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', alpha=0.1, elinewidth=3, 
                capsize=0, zorder=0, ms=5, mfc='blue', mec='blue', ecolor='blue');
    # Highlighting residue in question
    ax1.errorbar(a[0], b[0], xerr=a[1], yerr=b[1], fmt='o', ms=10, mfc='red', mew=0, 
                elinewidth=10, ecolor='red');
    ax1.set_ylabel(r"Major Population", fontsize=20);
    ax1.set_title(r"Res: %s, $k_{ex} vs. P_A$" %( res), fontsize=20);
    ymin, ymax = .45, 1.05
    ax1.set_ylim(ymin, ymax);

    if -50. < a[0] < 1500.:
        xmin, xmax = -50, 1500
        ax1.set_xlim(xmin, xmax);
        marker = [a[0], b[0]]
    elif np.isnan(a[0]) or np.isnan(b[0]):
        xmin, xmax = -50, 1500
        ax1.set_xlim(xmin, xmax);
        marker = [600, .6]
    else:
        xmin, xmax = -50, np.round(a[0] + 300, 0)
        ax1.set_xlim(xmin, xmax);
        marker = [a[0], b[0]]
    if hist == True:
        bins=15
        ax2.hist(x, bins=bins, range=(0, xmax))
        ax2.set_xlim(xmin, xmax)
        ax3.hist(y, bins=bins, range=(.5, 1), orientation='horizontal')
        ax3.set_ylim(ymin, ymax)
        ax2.set_xlabel(r"$k_{ex} (\frac{1}{s})$", fontsize=20);
    else:
        ax.set_xlabel(r"$k_{ex} (\frac{1}{s})$", fontsize=20);
    
    dw = df.loc[res,'dw_value':'dw_error']
    chi = df.loc[res, 'chi2_value']
    string = "$k_{ex,%i}$ = %.2f$\pm$%.2f\n$P_{A,%i}$ = %.2f$\pm$%.2f\n$\Delta\delta_{%i}$ = %.2f$\pm$%.2f\
    \n$\chi^2_{%i}$ = %.2f"
    ax1.annotate(string %(res, a[0], a[1], res, b[0], b[1], res, dw[0], dw[1], res, chi), 
                xy=(marker[0], marker[1]), xytext=(.3, .2),
                arrowprops=dict(facecolor='black', lw=2, arrowstyle='->'), 
                textcoords='axes fraction', fontsize=20, 
               )
    plt.tight_layout()
    return