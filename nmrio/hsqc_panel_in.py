import os
import matplotlib as mpl

def hsqc_panel_in(img_dir):
    """
    Description:  Given a directory containing only panels of an HSQC spectrum, preferably
    six panels, this function imports each image using matplotlib and returns a list of 
    ndarrays.
    
    Input:  img_dir; image directory
    """
    
    
    # Image import: HSQC Panels
    hsqc_imgs = list()
    for f in os.listdir(img_dir):
        if '.png' in f:
            hsqc_imgs.append(mpl.image.imread(os.path.join(img_dir, f)))
        
    return hsqc_imgs
