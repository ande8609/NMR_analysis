import matplotlib.pyplot as plt
import rd_curve, rd_residuals

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
    ax2.set_xlabel(r"$\nu_{CPMG} (Hz)$", fontsize=20);
    ax2.text(1100, 0.0, 'Residuals', transform=ax2.transData, fontsize=12)
    plt.xlim(-100, xmax);
    plt.tight_layout();
    return