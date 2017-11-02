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