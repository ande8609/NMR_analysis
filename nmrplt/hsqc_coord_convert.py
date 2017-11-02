from .hsqc_img_select import hsqc_img_select
import numpy as np

def hsqc_coord_convert(h, n, hsqc_imgs, hzoom=1, nzoom=1, hoff=0, noff=0):
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
    
    h_ex, n_ex = .55, 1
    h_avg, n_avg = np.mean(h)+hoff, np.mean(n)+noff
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
