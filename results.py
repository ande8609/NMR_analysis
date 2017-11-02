import pandas as pd
import os
from . import cpmg_plotting_fxns as cpf

#########################
# RELAX PAR/DATA IMPORT #
#########################
class results_relax(object):
    """
    A parent class to hold results from fits to RD data with relax.
    
    Subclasses will be for individual models with different parameter names.
    
    TO DO: 
    Nothin' yet.
    """
    
    
    
    def __init__(self, path):
        # Hiding raw data and parameters, they be nasty
        self._data, self._pars = self.load_data(path)
        # To keep track of where the data comes from 
        self.dpath = path
        self.parameters = list(self._pars.columns)
        self.residues = list(self._pars.index)
        # Data by field, useful for plotting?
        try:
            self.d600 = self._data.xs(600.53, level=1)
        except:
            self.d600 = None
        try:
            self.d800 = self._data.xs(799.70, level=1)
        except:
            self.d800 = None
        # Parameter setting, could use more foo (fu?)
        for par in self.parameters:
            value = self._pars[par]
            names = par.split("_")
            if len(names) > 2:
                name = "_".join([names[0], names[1], names[4].split('.')[0]])
            else:
                name = "_".join(names)
            setattr(self, name, value)
        # Calculate the Rex for each field and set attrs
        self._Rex_calc()
        
    def load_data(self, path):
        """
        Uses functions defined in the cpmg_plotting_fxns module of the rd_fxns package
        to load RD data and fit parameters from the given path.
        """
        # Loading raw parameters and data (including N/A values)
        par = cpf.read_param_dir(path)
        dat = cpf.read_relax_dir(path)
        # Getting rid of the N/A values
        par.dropna(axis=1, how='all', inplace=True)
        par.dropna(axis=0, thresh=5, inplace=True)
        dat.dropna(axis=1, how='all', inplace=True)
        return dat, par
    
    def by_res(self, res):
        """
        Pulls data from a big nasty multiindexed dataframe by residue. Returns two 
        series of data, the first collected at 600 MHz and the second at 800 MHz.
        """
        
        if res in self.residues:
            try:
                _600 = self.d600.xs(res, level=0)['R2eff_fit']
            except:
                _600 = None
            try:
                _800 = self.d800.xs(res, level=0)['R2eff_fit']
            except:
                _800 = None
            return _600, _800
        else:
            print("Residue not in data.")
            return
    
    def _Rex_calc(self):
        _600_list, _800_list = list(), list()
        
        for res in self.residues:
            try:
                _rex600 = self.by_res(res)[0].iloc[0] - self.r20_SQ_600[res]
                _rex800 = self.by_res(res)[1].iloc[0] - self.r20_SQ_799[res]
            except:
                _rex600 = None
                _rex800 = None
            _600_list.append(_rex600)
            _800_list.append(_rex800)
        self.rex600 = pd.Series(_600_list, index=self.residues)
        self.rex800 = pd.Series(_800_list, index=self.residues)




class results_relax_CR72(results_relax):
    """
    A class to hold results from a CR72 fit to RD data with relax.
    
    TO DO: 
    Figure out how to parse the multiindexed self._data pandas dataframe, its kind of a monster.
    """
    
    def __init__(self, path):
        results_relax.__init__(self, path)
        # Parsing parameters
        self.dw = self._pars[['dw_value', 'dw_error']]
        self.kex = self._pars[['kex_value', 'kex_error']]
        self.pa = self._pars[['pA_value', 'pA_error']]
        self.chisqr = self._pars['chi2_value']
        self.chi_reduced = self.chisqr/(len(self._data.index.get_level_values(2).unique()) - 5)
        # Parsing Data, 
        self.data = None
        self.fit = None
        
    def cull_data(self, chiCut=75, dwCut=.2, kexCut=1.):
        chi_set = [res for res, chi in self.chisqr.iteritems() if chi < chiCut]
        dw_set1 = [res for res in self.residues if self.dw['error'][res]/self.dw['value'][res] < dwCut]
        dw_set2 = [res for res in self.residues if self.dw['value'][res] > .35]
        kex_set = [res for res in self.residues if self.kex['value'][res] > kexCut]
        #dw_err_set = [res for res, err in self.dw['error'].iteritems() if err < dw_errCut]
        self.slice = sorted(set(chi_set) & set(dw_set1) & set(dw_set2) & set(kex_set))
        return self.slice

#####################
# LMFIT JSON IMPORT #
#####################
class results_lmfit(object):
    """
    A parent class to import and hold results in JSON format from fits to RD data 
    performed with the module lmfit.
    
    Subclasses will be for individual models with different parameter names. Done 
    so far are:  TSMFK01
    
    TO DO: 
    Build a subclass for CR72 data fit with lmfit.
    """
    
    def __init__(self, path):
        self._results_df = self.read_json_dir(path)
        self.dpath = path
        self.parameters = list(self._results_df.index.get_level_values(1).unique())
        self.residues = list(self._results_df.index.get_level_values(0).unique())
        self.slice = None
        try:
            self.chisqr = pd.Series.from_csv(os.path.join(path, 'chisqr.csv'))
        except:
            print("No chi squared .csv file found in dpath.")
    
    def load_data(self, new_path=None):
        if new_path != None:
            self.dpath = new_path
            self.read_json_dir(self.dpath)


    def read_json(self, fpath):
        """
        Read a TSMFK01 fit json exported by lmfit using pandas.
        """
        json_cols = ['par', 'value', 'hold', 'nan', 'lower', 
                     'upper', 'error', 'cross_corr', 'step?']
        test_df = pd.read_json(fpath)
        test_df.columns = json_cols
        test_df.index = test_df.par
        test_df.drop(['par', 'hold', 'nan', 'lower', 'upper', 'step?'], axis=1, inplace=True)
        return test_df
    
    def read_json_dir(self, dpath):
        """
        Read a directory full of json files exported by lmfit for TSMFK01
        fits to CPMG data.
        """
        results_list = list()
        resnums = list()
        for f in os.listdir(dpath): 
            if '.json' in f:
                resnums.append(int(f.split('.')[0].split('_')[1]))
                results_list.append(self.read_json(os.path.join(dpath, f)))
        results_df = pd.concat(results_list, keys=resnums, names=['resnum'])
        return results_df


class results_lmfit_TSMFK01(results_lmfit):
    
    def __init__(self, path):
        results_lmfit.__init__(self, path)
        self.dw = self._results_df.xs('dw_kay', level=1)
        self.kex = self._results_df.xs('kab', level=1)
        self.r20a = self._results_df.xs('r20a1', level=1)
        self.r20b = self._results_df.xs('r20a2', level=1)
        self.slice = None
        try:
            self.chisqr = pd.Series.from_csv(os.path.join(path, 'chisqr.csv'))
            self.chi_reduced = self.chisqr / (26-4)
        except:
            self.chisqr = None
            self.chi_reduced = None
            print("No chi squared .csv file found in dpath.")
      
    def cull_data(self, chiCut=75, dwCut=.2, kexCut=1.):
        chi_set = [res for res, chi in self.chisqr.iteritems() if chi < chiCut]
        dw_set1 = [res for res in self.residues if self.dw['error'][res]/self.dw['value'][res] < dwCut]
        dw_set2 = [res for res in self.residues if self.dw['value'][res] > .35]
        kex_set = [res for res in self.residues if self.kex['value'][res] > kexCut]
        #dw_err_set = [res for res, err in self.dw['error'].iteritems() if err < dw_errCut]
        self.slice = sorted(set(chi_set) & set(dw_set1) & set(dw_set2) & set(kex_set))
        return self.slice