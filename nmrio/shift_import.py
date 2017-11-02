
def shift_import(path, isdir=True):
    """
    Description:  A function that imports a peak list or lists and returns the 1H and 15N 
    chemical shift lists.
    
    Input:  path - path to the peak list or directory
    isdir; boolean True/False, true if a directory is passed, false if not.
    """
    
    if isdir == False:
        big_df = nmrio.ccp_tsv_in(path)
        h1_shift = big_df.xs('Position F1')
        n15_shift = big_df.xs('Position F2')
    else:
        big_df = nmrio.read_tsv_peak_dir(path)
        h1_shift = big_df.xs('Position F1', level='feature')
        n15_shift = big_df.xs('Position F2', level='feature')
    return h1_shift, n15_shift
