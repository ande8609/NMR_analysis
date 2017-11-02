from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import re
from scipy import stats
import fnmatch

"""
Note that this file is depricated as of 20160620.  Use the modules in the rd_fxns 
directory instead.
"""
##################################################################
##################################################################
##################################################################
# Functions from CPMG HSQC Zoom ##################################
##################################################################
##################################################################
##################################################################

def rex_csv(file):
    
    import pandas as pd
    df = pd.read_csv(file,
            header=None,names=['Rex'],index_col=0,skiprows=2).dropna().transpose()
    return df


def read_sparky(file):
    
    import pandas as pd
    #df = pd.read_table(file, header=None, names=['Ass','w1','w2','ht','vol'],delim_whitespace=True,na_values='?-?').dropna()
    df = pd.read_table(file, delim_whitespace=True, na_values='?-?').dropna()
    #df[['w1','w2','ht','vol']] = df[['w1','w2','ht','vol']].astype(float)
    #df = df.set_index('Ass').transpose()
    resnames = [int(re.findall(r'\d+',i)[0]) for i in df['Assignment']]
    df.index = resnames
    df = df.sort_index(axis=0)
    return df

def read_tsv_peak_list(f):
    
    import pandas as pd
    df = pd.read_csv(f,sep='\t',usecols=[2,3,5,6,7,8,9],index_col=2)
    index = [int(re.findall(r'\d+',i)[0]) for i in df.index.values]
    AA = [re.findall(r'[A-Za-z]{1,3}',i)[0] for i in df.index.values]
    df['AA'] = pd.Series(AA)
    df.index = index
    return df.transpose().sort_index(axis=1)

"""
def shift_import(path, isdir=True):
    if isdir == False:
        big_df = read_tsv_peak_list(path)
        h1_shift = big_df.xs('Position F1')
        n15_shift = big_df.xs('Position F2')
    else:
        big_df = read_tsv_peak_dir(path)
        h1_shift = big_df.xs('Position F1', level='feature')
        n15_shift = big_df.xs('Position F2', level='feature')
    return h1_shift, n15_shift
"""

def shift_import(path, isdir=True):
    """
    A function to import peak lists and return 1H and 15N 
    chemical shift lists.
    """
    
    if isdir == False:
        df = read_sparky(path)
        h1_shift = df['w1']
        n15_shift = df['w2']
    else:
        big_df = read_sparky(path)
        h1_shift = big_df.xs('w1', level='feature')
        n15_shift = big_df.xs('w2', level='feature')
    return h1_shift, n15_shift


def read_tsv_peak_dir(data_dir):
    """
    A function to read all peak list files of the form '*list.tsv' in the 
    given data directory.  This version also has functionality hard coded
    in it that pulls the spectrometer frequency out, which is predicated
    on the filename starting with that frequency.
    """
    import os
    import fnmatch
    import re
    
    # Finding all peak list files, *list.tsv, in the current directory
    # Also pulling out the 
    df_list = list()
    keys = list()
    for f in os.listdir(data_dir):
        if fnmatch.fnmatch(f, '*list.tsv'):
            df_list.append(read_tsv_peak_list(os.path.join(data_dir, f)))
            keys.append(re.findall(r'\d+', f)[0])
    
    # Making a large dataframe from all peak list files
    big_df = pd.concat(df_list,axis=0,keys=keys,names=['peak list','feature'])
    return big_df

def img_import(img_dir):
    """
    
    """
    
    import os
    import matplotlib as mpl

    # Image import: HSQC Panels
    hsqc_imgs = list()
    for f in os.listdir(img_dir):
        hsqc_imgs.append(mpl.image.imread(os.path.join(img_dir, f)))
        
    return hsqc_imgs

def hsqc_img_select(h, n):
    """
    Given a proton and nitrogen chemical shift pair, this function chooses
    the image corresponding to that peak.
    
    Functionality to choose two image panels and may be usefule and could
    be added in later by modifying the conditionals and making a parser
    function elsewhere.
    """
    if n < 118.75:
        if h > 8.67:
            img = 0
        if h > 7.33 and h < 8.67:
            img = 1
        if h < 7.33:
            img = 2
    if n > 118.75:
        if h > 8.67:
            img = 3
        if h > 7.33 and h < 8.67:
            img = 4
        if h < 7.33:
            img = 5
    return img

def hsqc_coord_convert(h, n, hsqc_imgs, hzoom=1, nzoom=1):
    """
    This function takes a pair of proton, nitrogen chemical shift frequencies
    and calculates the conversion from chemical shift units to pixels from 
    the edges of the image panel to that peak.
    
    Inputs:
    h, n:  single values or vectors representing proton and nitrogen chemical
    shifts of a single or multiple peaks.
    hzoom, nzoom:  The factor to zoom in, or out, beyond the set amount so that
    zooming can be done on the fly.
    
    Outputs are:
    xmax, xmin:  Max & min in pixels, for image cropping, use ax.set_xlim()
    ymax, ymin:  Max & min in pixels, for image cropping, use ax.set_ylim()
    hmax, hmin:  Max & min in 1H freq., for image tick labels
    nmax, nmin:  Max & min in 15N freq., for image tick labels
    img_num:     Panel chosen with hsqc_img_select
    
    Note that h_ex and n_ex are factors to choose the degree of zoom in each
    dimension.  Averages are taken 
    """
    
    import numpy as np
    h_ex, n_ex = .55, 1
    h_avg, n_avg = np.mean(h), np.mean(n)
    img_num = hsqc_img_select(h_avg,n_avg)
    if img_num == 0:
        hmax, hmin, nmax, nmin = 10.0, 8.67, 118.75, 104.5
    if img_num == 1:
        hmax, hmin, nmax, nmin = 8.67, 7.33, 118.75, 104.5
    if img_num == 2:
        hmax, hmin, nmax, nmin = 7.33, 6.0, 118.75, 104.5
    if img_num == 3:
        hmax, hmin, nmax, nmin = 10.0, 8.67, 133.0, 118.75
    if img_num == 4:
        hmax, hmin, nmax, nmin = 8.67, 7.33, 133.0, 118.75
    if img_num == 5:
        hmax, hmin, nmax, nmin = 7.33, 6.0, 133.0, 118.75
    # Conversions in pixels/ppm
    hconvert, nconvert = hsqc_imgs[img_num].shape[1] / (hmax - hmin), hsqc_imgs[img_num].shape[0] / (nmax - nmin)
    
    x = (hmax - h_avg) * hconvert
    y = (n_avg - nmin) * nconvert
    
    xmax = int(np.around((hmax - np.max(h) + h_ex/hzoom) * hconvert, 0))
    xmin = int(np.around((hmax - np.min(h) - h_ex/hzoom) * hconvert, 0))
    ymin = int(np.around((np.max(n) - nmin - n_ex/nzoom) * nconvert, 0))
    ymax = int(np.around((np.min(n) - nmin + n_ex/nzoom) * nconvert, 0))
    hmin = (np.max(h) + h_ex)
    hmax = (np.min(h) - h_ex)
    nmin = (np.max(n) - n_ex)
    nmax = (np.min(n) + n_ex)
    return [xmax, xmin, ymax, ymin, hmax, hmin, nmax, nmin, img_num, x, y]

def hsqc_zoom_subplot(res=None, ax=None, hsqc_imgs=None, h1_shift_df=None, n15_shift_df=None, hzoom=None, nzoom=None):
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
    """
    
    import matplotlib.pyplot as plt
    if res == None:
        res = 152
    if ax == None:
        fig, ax = plt.subplots(1, 1, figsize = (10, 7.5));
    if hzoom == None or nzoom == None:
        hzoom, nzoom = 2.25, .75
    h = h1_shift_df[res]
    n = n15_shift_df[res]
    xmax, xmin, ymax, ymin, hmax, hmin, nmax, nmin, img_num, h_mean, n_mean = hsqc_coord_convert(h,n, hsqc_imgs,hzoom, nzoom)
    ax.imshow(hsqc_imgs[img_num], zorder=0);
    ax.plot(np.repeat(h_mean, 10), np.linspace(ymax, ymin, 10), '--', color='grey', lw=5, zorder=1);
    ax.plot(np.linspace(xmin, xmax, 10), np.repeat(n_mean, 10), '--', color='grey', lw=5, zorder=1);
    #ax1.scatter(h_mean, n_mean, color='grey', s=150);
    ax.add_patch( mpl.patches.Rectangle( (h_mean+(37.5), n_mean-(100)), 165, 75, 
                                         fill=False, ec='grey', lw=5) )
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


################################################################################
################################################################################
################################################################################
# Functions from CPMG Data Subplotting #########################################
################################################################################
################################################################################
################################################################################

# Data import functions

def read_relax_file(directory, filename):
    """
    A function that takes a directory and a 'disp*.out' file from a relax analysis and returns
    a dataframe multi-indexed by field strength at the first level and nu_cpmg at the second
    level.
    """
    df = pd.read_table( os.path.join(directory, filename), skiprows=1, usecols=range(2,7), header=None, 
                             names=['Field','nu_cpmg','R2eff_fit','R2eff_calc','R2eff_error'], 
                             delim_whitespace=True, na_values='-', index_col=[0,1] )
    return df

def read_relax_dir(directory):
    """
    Given a relax analysis directory containing 'disp*.out' files, this function reads all such
    files and returns a multi-indexed dataframe with the first level as residue number ('Residue') 
    the second level field strength ('Field') and the third index level cpmg pulsing frequency
    ('nu_cpmg').
    """
    import os
    import fnmatch
    import re
    
    res_list = list()
    df_list = list()
    for f in os.listdir(directory):
        if fnmatch.fnmatch(f, "disp*.out"):
            res_list.append(int(re.findall("[0-9]+",f)[0]))
            df_list.append( read_relax_file(directory, f) )
    
    big_df = pd.concat( df_list, axis=0, keys=res_list, names=['Residue', 'Field', 'nu_cpmg'] ).sort_index()
    # Needed for quick retrieval when slicing
    big_df.sortlevel(inplace=True)
    return big_df

# Data Plotting Functions

def rd_curve_plot(df, res=None, field=None, ax=None, typ=None):
    """
    This function plots the RD data and fit points.
    
    INPUTS:
    df:  multi-indexed dataframe
    res:  residue number
    field: 0 or 1
    ax:  subplotting axis to be used
    typ:  'data' or 'fit'
    """
    
    import matplotlib.pyplot as plt
    import pandas as pd
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

def rd_residuals_plot(df, res=None, field=None, ax=None):
    """
    Plots the residuals of the fit (y(data) - y(fit)).
    
    INPUTS:
    df:  multi-indexed dataframe
    res:  residue number
    field:  0 or 1, corresponding to multi-index level 1
    ax:  subplotting axis to be used
    """
    
    import matplotlib.pyplot as plt
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

def rd_agg_plot(df, res=None, ax1=None, ax2=None ):
    """
    Plots the RD curves and residuals together.
    
    INPUTS:
    df:  multi-indexed dataframe
    res:  residue number
    axs:  two subplotting axes to be used, one for curves
          one for residuals.
    """
    
    import matplotlib.pyplot as plt
    if res == None:
        res = 152

    if ax1 == None or ax2 == None:
        ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=3);
        ax2 = plt.subplot2grid((4, 1), (3, 0), sharex=ax1);
    
    xmax = 1300
    # First Axes: Data + Fit
    rd_curve_plot(df, res, typ='data', ax=ax1, field=1)
    rd_curve_plot(df, res, typ='fit', ax=ax1, field=1)
    rd_curve_plot(df, res, typ='data', ax=ax1, field=0)
    rd_curve_plot(df, res, typ='fit', ax=ax1, field=0)
    ax1.legend(framealpha = 1., loc=7);
    ax1.set_ylabel(r"$R_{2,eff} (\frac{1}{s})$", fontsize=20)
    ax1.set_title("Residue %s, RD Curves" %(res), fontsize=20)
    
    # Second Axes: Residuals
    rd_residuals_plot(df, res, ax=ax2, field=0)
    rd_residuals_plot(df, res, ax=ax2, field=1)
    ax2.set_xlabel(r"$\tau_{CPMG} (ms)$", fontsize=20);
    ax2.text(1100, 0.0, 'Residuals', transform=ax2.transData, fontsize=12)
    plt.xlim(-100, xmax);
    plt.tight_layout();
    return

################################################################################
################################################################################
################################################################################
# Functions from CPMG Simulate over Data #######################################
################################################################################
################################################################################
################################################################################


def cr72_simulate(r20a=20,r20b=20,pA=.95,kex=500.,dw_h=.1,dw_n=1.):
    """
    DESCRIPTION:
    Function that uses the `relax` NMR RD modules to calculate the RD curve
    for a given set of parameters for both 600 and 800 MHz and then plots 
    both curves.
    
    INPUTS:
    r20a, r20b:  intrinsic r2 rates for states a and b, respectively, usually
                 the same.
    pA:  population of the major state.
    kex:  k1 + k2, the sum of the forward and reverse rates of exchange.
    dw_h, dw_n:  differences in chem. shift between the two states in ppm for 
    1H and 15N, respectively
    """
    
    
    import sys
    sys.path.append('/usr/local/relax-4.0.2/lib')
    sys.path.append('/usr/local/relax-4.0.2')
    import relax
    import nmr
    from dispersion.cr72 import r2eff_CR72
    from ipywidgets import interact, interactive, fixed
    import ipywidgets as widgets
    import matplotlib.pyplot as plt
    import numpy as np

    # Choosing the major population slider
    if pA == None:
        pA = 0.95
    else:
        pA = pA
        
    # Setting the intrinsic r2 value, held constant
    if r20a == None or r20b == None:
        r2int=20.
    elif r20b == None:
        r2int=r20a
    else:
        r2int=r20b
    
    # Setting kex
    if kex == None:
        kex = 500.
    else:
        kex = kex
        
    # Setting dw_h & dw_n
    if dw_h == None or dw_n == None:
        dw_h = .1
        dw_n = 1.
    else:
        dw_h = dw_h
        dw_n = dw_n

    fields = [600e6, 800e6]
    colors = ['blue', 'red']
    labels = ['600 MHz', '800 MHz']
    r20a_orig = r20a
    r20b_orig = r20b
    cpmg_frqs = np.linspace(20,1200,150)
    fig, ax = plt.subplots(1, 1, figsize=(12,10))
    
    for n, field in enumerate(fields):
        # Convert dw values from ppm to rad/s for both 1H & 15N
        dw1_con = nmr.frequency_to_rad_per_s(dw_h, field, '1H')
        dw2_con = nmr.frequency_to_rad_per_s(dw_n, field, '15N')
        dw = dworig = ((dw1_con)**2 + (dw2_con/6.51)**2)**(.5)
        dw_ppm = ((dw_h)**2 + (dw_n/6.51)**2)**(.5)
        back_calc = np.zeros(cpmg_frqs.shape)            
        back_calc = np.zeros(cpmg_frqs.shape)
        k = kex*2.*np.pi
        # Use the native `relax` CR72 function, stores values in back_calc
        r2eff_CR72(r2int, r2int, r2int, r2int, pA, dw, dworig, k, cpmg_frqs, back_calc)
        ax.plot(cpmg_frqs, back_calc, lw=10, alpha=.75, color=colors[n])
    print('kex=%i, dw=%.2f' %(int(kex), dw_ppm))
    plt.legend(labels);
    plt.ylabel("R$_2$ s$^{-1}$");
    plt.xlabel("CPMG Frequency Hz")
    return

def cr72_manipulate(pA=None, r20a=None, r20b=None, kex=None, dw_h=None, dw_n=None, num=None):
    """
    My attempt at using interact in a function.
    
    INPUTS:
    r20a, r20b:  intrinsic r2 rates for states a and b, respectively, usually
                 the same.
    pA:  population of the major state.
    kex:  k1 + k2, the sum of the forward and reverse rates of exchange.
    dw_h, dw_n:  differences in chem. shift between the two states in ppm for 
    1H and 15N, respectively


    
    TO DO:
    Add in the ability to import parameter values from the `relax` analysis
    for a particular residue and set the values and ranges in the interact 
    call automatically.
    
    Make another function to plot RD data for a given residue and then a third
    function that calls both the data and simulation functions.  I MAY HAVE THE 
    FIRST FUNCTION HERE IN MY RD PLOTTING FROM BEFORE.
    
    Add in functionality to calculate whether a residue is in slow, fast, or 
    intermediate exchange.  That could be a function itself.
    """

    
    import sys
    sys.path.append('/usr/local/relax-4.0.2/lib')
    sys.path.append('/usr/local/relax-4.0.2')
    import relax
    import nmr
    from dispersion.cr72 import r2eff_CR72
    from ipywidgets import interact, interactive, fixed
    import ipywidgets as widgets
    import matplotlib.pyplot as plt
    import numpy as np

    
    # The number of steps for each slider
    if num == None:
        num = 15
    
    # Choosing the major population slider
    if pA == None:
        lower, upper = .5, 1.
        step, val = (upper-lower)/num, 0.9
        pA_slider=widgets.FloatSlider(min=lower, max=upper, step=step, value=val)
    else:
        lower, upper = .5, 1.
        step, val = (upper-lower)/num, pA
        pA_slider=widgets.FloatSlider(min=lower, max=upper, step=step, value=val)
    
    # Setting the intrinsic r2 value, held constant
    if r20a == None or r20b == None:
        r2int=fixed(20)
    elif r20b == None:
        r2int=fixed(r20a)
    else:
        r2int=fixed(r20b)
    
    # Setting kex
    if kex == None:
        lower, upper, step, val = 100., 600., (600.-100.)/num, 500.
        kex_slider = widgets.FloatSlider(min=lower, max=upper, step=step, value=val)
    else:
        lower, upper = .1 * kex, 1.5 * kex
        step, val = (upper-lower)/num, kex
        kex_slider = widgets.FloatSlider(min=lower, max=upper, step=step, value=val)
        
    # Setting dw_h & dw_n
    if dw_h == None or dw_n == None:
        lower, upper = .1, 4.
        step, val = (upper-lower)/num, 1.
        dw_n_slid = widgets.FloatSlider(min=lower, max=upper, step=step, value=val)
        lower, upper = .01, 1
        step, val = (upper-lower)/num, .1
        dw_h_slid = widgets.FloatSlider(min=lower, max=upper, step=step, value=val)

    else:
        lower, upper = dw_h * .1, dw_h * 1.5
        step, val = (upper-lower)/num, dw_h
        dw_h_slid = widgets.FloatSlider(min=lower, max=upper, step=step, value=val)
        lower, upper = dw_n * .1, dw_n * 1.5
        step, val = (upper-lower)/num, dw_n
        dw_n_slid = widgets.FloatSlider(min=lower, max=upper, step=step, value=val)
    
    
    interact(cr72_simulate, 
             pA=pA_slider, 
             r20a=r2int, r20b=r2int, 
             kex=kex_slider, 
             dw_h=dw_h_slid, 
             dw_n=dw_n_slid);
    return

################################################################################
################################################################################
################################################################################
# Functions from Fit Parameter Plotting ########################################
################################################################################
################################################################################
################################################################################

# Data import functions

def read_param_dir(data_dir):
    """
    This was a pain in the ass.
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

# Data Plotting Functions

def kex_pa_plot(df, res=None, ax1=None, ax2=None, ax3=None, hist=True):
    """
    Plotting kex vs. pA CPMG parameters.
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


def img_read(path):
    """
    Reads images generated by CPMG plotting functions.
    INPUTS:
    path:  The relative path to the directory with all the
           the images.
    """
    import progressbar as pb
    
    names = list()
    objects = list()
    flist = os.listdir(path)
    
    # Initializing the progress bar
    num = len(flist)
    widgets = ["    ", pb.Percentage(), "    ", pb.Bar(), "    ", pb.Timer(), "    ", pb.ETA(), "    "]
    # Define bar
    bar = pb.ProgressBar(maxval=num, widgets=widgets)
    # Start bar
    bar.start()
    
    for i, f in enumerate(flist):
        if 'png' in f:
            names.append(int(re.findall("[0-9]+",f)[0]))
            objects.append(plt.imread(os.path.join(path, f)))
            # Update bar
            bar.update(i+1)
    # Finish bar
    bar.finish()
    ser = pd.Series(objects, index=names)
    return ser

def plot_figures(figures, nrows=2, ncols=2, save=False):
    """
    
    """
    
    fig, axeslist = plt.subplots(ncols=ncols, nrows=nrows, figsize = (11, 8.5))
    for col in range(len(figures)):
        axeslist.ravel()[col].imshow(figures[col])
        axeslist.ravel()[col].set_axis_off()
    fig.subplots_adjust(hspace=0, wspace=0, top=1, bottom=0, left=0, right=1)

def save_figs(path, df, seq):
    """
    
    """
    
    import progressbar as pb
    import os
    cols = df.columns

    # Initializing the progress bar
    num = len(cols)
    widgets = ["    ", pb.Percentage(), "    ", pb.Bar(), "    ", pb.Timer(), "    ", pb.ETA(), "    "]
    # Define bar
    bar = pb.ProgressBar(maxval=num, widgets=widgets)
    # Start bar
    bar.start()

    
    for i, col in enumerate(cols):
        plot_figures(df[col], nrows=2, ncols=2);
        plt.title('%s%s' %(col, seq[i]), fontsize=20)
        plt.savefig(os.path.join(path, '%s.png' %(col)), 
                    type='png', dpi=200)
        plt.close()
        bar.update(i+1)
    bar.finish()

def img_df(hsqc_imgs, curve_imgs, param_imgs, struc_imgs):
    sindex = sorted(list(set(curve_imgs.index.values) & 
                  set(param_imgs.index.values) & 
                  set(hsqc_imgs.index)))
    index = ['curve', 'params', 'hsqc', 'struc']
    df = pd.DataFrame([curve_imgs[sindex], param_imgs[sindex], 
                       hsqc_imgs[sindex], struc_imgs[sindex]], 
                       index=index)
    return df

def img_df_unsigned(hsqc_imgs, curve_imgs, param_imgs):
    sindex = sorted(list(set(curve_imgs.index.values) & 
                  set(param_imgs.index.values) & 
                  set(hsqc_imgs.index)))
    index = ['curve', 'params', 'hsqc']
    df = pd.DataFrame([curve_imgs[sindex], param_imgs[sindex], 
                       hsqc_imgs[sindex]], 
                       index=index)
    return df