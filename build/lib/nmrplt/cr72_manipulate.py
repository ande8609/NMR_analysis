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
