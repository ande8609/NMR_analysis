import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

def slev1_big(big_df, feature='Height', level=1):
    """
    A function for slicing the monsters I've created.  Takes in a 
    large multi-indexed dataframe with resnums as columns and a 
    level zero index labeled by peak list, temps or rd freq, and 
    a level one index labeled by feature, which is what I'm slicing
    here.  Returns a dataframe sliced according to some level 1 
    feature.
    """
    
    level_names = list(set(big_df.index.get_level_values(level)))
    try:
        df = big_df.xs(feature, level=level)
        return df
    except:
        print('Error:  feature passed is not in the level, choose one of the following:')
        print(level_names)
    return

def int_changes(res, big_df, volume=True, ax=None, slow=False):
    """
    A function to plot the peak intensity changes as a function
    of some variable, usually temp. or binding partner.
    """
    
    size = 250
    if ax == None:
        ax = plt.subplot2grid((1, 1), (0, 0));
    temps = sorted(list(set(big_df.index.get_level_values(0))))[0::2]
    x = [int(re.findall(r'\d+', i)[0]) for i in temps]
    # Getting res name
    resn = list(set(slev1_big(big_df, 'resn')[res]))[-1]
    # Volume vs. height
    if not volume:
        y = slev1_big(big_df)
    else:
        y = slev1_big(big_df, 'Volume')
    # Choosing folded vs. unfolded for V135G   
    names = ['Folded', 'Unfolded']
    colors = ['blue', 'red']
    if slow:
        num = 2
    else:
        num = 1
    for i in range(num):
        y1 = y[res][i::2].values
        color = 'blue'
        ax.plot(x, y1, '.', c=colors[i], ms=25, label=names[i]);
        ax.plot(x, y1, '--', zorder=0, lw=5, c=colors[i], alpha=.5)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.set_xlim(9, 26);
    ax.legend(loc='lower center', framealpha=.5)
    ax.set_xlabel('Temperature ($\degree$C)')
    ax.set_ylabel('Peak Intensity')
    ax.set_title('%s%s $\Delta$Intensity' %(resn, res))
    plt.tight_layout()


def lw_changes(res, big_df, ax1=None, ax2=None, slow=False):
    """
    A function to plot the peak intensity changes as a function
    of some variable, usually temp. or binding partner.
    """
    
    size = 250
    if ax1 == None and slow == True:
        #fig = plt.figure(figsize=(10, 8))
        ax1 = plt.subplot2grid((2, 1), (0, 0));
        ax2 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1);
        axs = [ax1, ax2]
    elif ax1 == None:
        #fig = plt.figure(figsize=(10, 8))
        ax1 = plt.subplot2grid((2, 1), (0, 0));
        axs = [ax1]
    temps = sorted(list(set(big_df.index.get_level_values(0))))[0::2]
    x = [int(re.findall(r'\d+', i)[0]) for i in temps]
    # Getting res name, lots of trouble to avoid nan
    resn = list(set(slev1_big(big_df, 'resn')[res]))[-1]
    # Getting both linewidths
    # Choosing folded vs. unfolded for V135G    
    names = ['Folded', 'Unfolded']
    
    if slow:
        num = 2
    else:
        num = 1
    for i in range(num):
        y1 = slev1_big(big_df, 'Line Width F1 (Hz)')[res][i::2]
        y2 = slev1_big(big_df, 'Line Width F2 (Hz)')[res][i::2]
        #y1 = y1
        color = 'red'
        axs[i].scatter(x, y1, c=color, s=size, label="$^1$H");
        axs[i].plot(x, y1, zorder=0, lw=5, c=color, alpha=.5, label='_')
        #y2 = y2.dropna(axis=1).iloc[1::2, :][res].values
        color = 'blue'
        axs[i].scatter(x, y2, c=color, s=size, label="$^{15}$N");
        axs[i].plot(x, y2, zorder=0, lw=5, c=color, alpha=.5, label='_')
        axs[i].set_ylabel('Line Width (Hz)')
        axs[i].set_title('%s%s $\Delta$Linewidth, %s' %(resn, res, names[i]))
        axs[i].legend(loc='lower right', framealpha=.5)
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    axs[i].set_xlabel('Temperature ($\degree$C)')
    plt.tight_layout()


def hsqc_changes(res, big_df, folded=True, ax1=None, ax2=None):
    """
    A function to plot the chemical shift changes as a function
    of some variable, usually temp. or binding partner.
    """
    
    if ax1 == None:
        #fig = plt.figure(figsize=(10, 8))
        ax1 = plt.subplot2grid((2, 1), (0, 0));
        ax2 = plt.subplot2grid((2, 1), (1, 0));
    resn = set(slev1_big(big_df, 'resn')[res]).pop()
    temps = sorted(list(set(big_df.index.get_level_values(0))))[0::2]
    temps = [int(re.findall(r'\d+', i)[0]) for i in temps]
    h1_shift_df = slev1_big(big_df, 'Position F1')
    n15_shift_df = slev1_big(big_df, 'Position F2')
    colors = ['blue', 'skyblue', 'plum', 'coral', 'red']
    #colors = ['purple', 'blue', 'green', 'orange', 'red']
    axs = [ax1, ax2]
    names = ['Folded', 'Unfolded']
    for i in range(2):
        x = h1_shift_df.dropna(axis=1).iloc[i::2, :][res].values
        y = n15_shift_df.dropna(axis=1).iloc[i::2, :][res].values
        for n, pair in enumerate(zip(x, y)):
            scat = axs[i].scatter(pair[0], pair[1], c=colors[n], 
                                  s=250, 
                                  label=str(temps[n])+'$\degree$C');
        axs[i].set_title('%s%s $\Delta\delta$, %s' %(resn, res, names[i]))
    ax1.invert_xaxis()
    ax1.invert_yaxis()
    ax2.invert_xaxis()
    ax2.invert_yaxis()
    plt.ticklabel_format(style='plain', axis='y')
    ax1.set_xlabel("$^1$H ppm")
    ax1.set_ylabel("$^{15}$N ppm")
    ax2.set_xlabel("$^1$H ppm")
    ax2.set_ylabel("$^{15}$N ppm")
    #fig.subplots_adjust(right=0.8)
    plt.legend()
    plt.tight_layout()


def vantHoff_df(big_df, volume=False):
    """
    Takes a multiindexed dataframe generated by ccp_tsv_in and 
    calculates Keq and ln(Keq) for slow exchange data.
    """
    
    if volume == True:
        intensity = 'Volume'
    else:
        intensity = 'Height'
    df_list = list()
    names = ['eq_uf', 'ln_uf', 'eq_fu', 'ln_fu']
    df_f_ht = slev1_big(big_df, intensity).iloc[::2,:-1]
    df_f_ht.index = [int(re.findall(r'\d+', i)[0]) for i in df_f_ht.index]
    df_u_ht = slev1_big(big_df, intensity).iloc[1::2,:-1]
    df_u_ht.index = [int(re.findall(r'\d+', i)[0]) for i in df_u_ht.index]
    df_list.append(df_u_ht.divide(df_f_ht).dropna(axis=1))
    df_list.append(df_list[0].applymap(np.log))
    df_list.append(df_f_ht.divide(df_u_ht).dropna(axis=1))
    df_list.append(df_list[2].applymap(np.log))
    big_df = pd.concat(df_list, keys=names, names=['Type', 'Temp'])
    return big_df.swaplevel(0, 1)


def vantHoff(res, big_df, uf=True, ln=True, ax=None, color=False):
    """
    A function that takes in a dataframe of ln(keq) values
    and a residue number and makes a van't Hoff plot in 
    addition to calculating the slope, which is -dH/R.
    """
    
    from scipy import stats
    if ax == None:
        fig, ax = plt.subplots(1, 1);
    
    big_vh = vantHoff_df(big_df)
    if ln == True and uf == True:
        df = slev1_big(vantHoff_df(big_df), 'ln_uf')
    elif ln == True and uf == False:
        df = slev1_big(vantHoff_df(big_df), 'ln_fu')
    elif ln == False and uf == True:
        df = slev1_big(vantHoff_df(big_df), 'eq_uf')
    elif ln == False and fu == True:
        df = slev1_big(vantHoff_df(big_df), 'eq_fu')
    else:
        print("Error, invalid slicing option")
    
    resn = list(set(slev1_big(big_df, 'resn')[res]))[-1]
    R = 1.9872035e-3 # in Kcal/mol
    if ln == True:
        xlab = r"$\frac{1}{T}$ $K^{-1}$"
        ylab = "ln(K$_{EQ}$)"
        x = 1./(df[res].index.values+273.3)
    else:
        xlab = "T ($\degree$C)"
        ylab = "K$_{EQ}$"
        x = (df[res].index.values)
    y = df[res].values
    m, b, r, p, stderr = stats.linregress(x, y)
    if color == True:
        temps = df[res].index.values
        colors = ['#0015e5', '#6500e9', '#e400ed', '#f1007a', '#f50500']
        for n, pair in enumerate(zip(x, y)):
            ax.scatter(pair[0], pair[1], s=250, c=colors[n], 
                       label=str(temps[n])+'$\degree$C')
    elif color == False:
        ax.scatter(x, y, s=250, label="Data")
    xtend = np.linspace(min(x)*.99, max(x)*1.01, 10)
    ax.plot(xtend, m*xtend + b, lw=5, zorder=0, alpha=.5, label="Fit", 
            c='black')
    ax.set_xlim(min(xtend), max(xtend))
    if ln == True:
        ax.text(.1,.1,"$\Delta$H=%.2f Kcal/mol\nr$^2$=%.3f\nstderr=%.2f" 
                %(-m*R, r**2, stderr),  
                transform=ax.transAxes, fontsize=20);
    plt.legend();
    ax.set_xlabel(xlab);
    ax.set_ylabel(ylab);
    ax.set_title("%s%s van't Hoff" %(resn, res))