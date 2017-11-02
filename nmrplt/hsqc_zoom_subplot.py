from .hsqc_coord_convert import hsqc_coord_convert
from .hsqc_img_select import hsqc_img_select
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def hsqc_zoom_subplot(res=None, hsqc_imgs=None, h1_shift_df=None, n15_shift_df=None, ax=None, hzoom=None, nzoom=None, hoff=0, noff=0):
    """
    Plots a region of an HSQC zoomed in on the residue in question.
    
    INPUTS:
    hsqc_imgs:  a list of six panels of the same hsqc spectra.  See hsqc_img_select
    and hsqc_coord_convert for panel numbering and coordinates.
    h1_shift_df:  dataframe containing 1H chemical shifts with columns
    labeled by amino acid number and rows as spec. freq. or temp. or anything.
    n15_shift_df:  dataframe containing 15N chemical shifts with columns
    labeled by amino acid number and rows as spec. freq. or temp. or anything.
    res:  residue number as an integer.
    ax:  axes object for subplotting.
    
    NOTES:  This should be changed to take in the big multiindexed data
    frame instead.
    """
    
    if res == None:
        res = 152
    if ax == None:
        fig, ax = plt.subplots(1, 1, figsize = (10, 7.5));
    if hzoom == None or nzoom == None:
        hzoom, nzoom = 2.25, .75
    h = h1_shift_df[res]
    n = n15_shift_df[res]
    xmax, xmin, ymax, ymin, hmax, hmin, nmax, nmin, img_num, h_mean, n_mean = hsqc_coord_convert(h,n, hsqc_imgs,hzoom, nzoom, hoff, noff)
    ax.imshow(hsqc_imgs[img_num], zorder=0);
    ax.plot(np.repeat(h_mean, 10), np.linspace(ymax, ymin, 10), '--', color='grey', lw=5, zorder=1);
    ax.plot(np.linspace(xmin, xmax, 10), np.repeat(n_mean, 10), '--', color='grey', lw=5, zorder=1);
    #ax1.scatter(h_mean, n_mean, color='grey', s=150);
    #ax.add_patch( mpl.patches.Rectangle( (h_mean+(37.5), n_mean-(100)), 165, 75, 
    #                                     fill=False, ec='grey', lw=5) )
    ax.set_xlim(xmin, xmax);
    ax.set_ylim(ymax, ymin);

    xticknum = len(ax.get_xticklabels());
    yticknum = len(ax.get_yticklabels());
    ax.set_xticklabels( np.around( np.linspace(hmin,hmax,xticknum), 1 ), fontsize = 20 );
    ax.set_yticklabels( np.around( np.linspace(nmin,nmax,yticknum), 1 ), fontsize = 20 );
    ax.set_xlabel('$^1$H ppm', fontsize=20);
    ax.set_ylabel('$^{15}$N ppm', fontsize=20);
    ax.text(.6,.2,"- 600 MHz",fontsize=20, backgroundcolor='white', transform=ax.transAxes);
    ax.text(.6,.12,"- 800 MHz",fontsize=20, backgroundcolor='white', transform=ax.transAxes);
    ax.text(.56,.2,"x",fontsize=36, backgroundcolor='white', transform=ax.transAxes);
    ax.text(.56,.12,"o",fontsize=36, color='red', backgroundcolor='white', transform=ax.transAxes);
    ax.set_title("Res. %i" %(res), fontsize = 20);
    plt.tight_layout()
    return