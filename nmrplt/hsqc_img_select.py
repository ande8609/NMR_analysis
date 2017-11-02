
def hsqc_img_select(h, n):
    """
    Given a proton and nitrogen chemical shift pair, this function chooses
    the image corresponding to that peak.
    
    Functionality to choose two image panels and may be usefule and could
    be added in later by modifying the conditionals and making a parser
    function elsewhere.
    """
    if n < 118.74:
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
