"""bez_models.py: Module is used to implement various sea-earth models and calculate TFs"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
os.system("rm -rf `find -type d -name .ipynb_checkpoints`:")

# Import required packages
import bezpy
import numpy as np
from scipy import constants as C
import pandas as pd

from multiprocessing import Pool
from fastkml import kml
from netCDF4 import Dataset

import plotlib

class OceanModel(object):
    """
    This class is to emulate electric and magnetic fields
    in a 1D ocean. We provide the 1D Earth conductivity model and
    one layer of ocean parameters (depth and resistivity) to emulate 
    1D ocean.
    """
    
    def __init__(self, model_name="BME", ocean_model={"depth":5e3, "rho":0.25}, flim=[1e-6, 1e0]):
        self.model_name = model_name
        self.ocean_model = ocean_model
        self.site = bezpy.mt.read_1d_usgs_profile("data/ocean_model_%s.txt"%model_name)
        self.freqs = np.linspace(flim[0], flim[1], int(flim[1]/flim[0])+1)
        self.functions = {
            "Ef2Bs": lambda Z, Zd, kd: Zd/(np.cosh(kd) + (Zd*np.sinh(kd)/Z)),
            "Bf2Bs": lambda Z, Zd, kd: 1./(np.cosh(kd) + (Zd*np.sinh(kd)/Z))
        }
        return
    
    def calcZ(self, layer="ocean"):
        if not hasattr(self, "Z"):
            freqs = np.copy(self.freqs)
            omega = 2*C.pi*self.freqs
            sigma_s = 1/self.ocean_model["rho"]
            k2 = 1.j*omega*C.mu_0*sigma_s
            k = np.sqrt(k2)
            Zo = (1.j*omega*C.mu_0/k)/(C.mu_0/1.e-3)
            Kf = self.site.calcZ(freqs)[1]
            self.Z = {"ocean": Zo, "floor": Kf}
        return self.Z[layer]
    
    def calcTF(self, kinds=["Ef2Bs", "Bf2Bs"]):
        omega = 2*C.pi*self.freqs
        Zd = self.calcZ("floor")
        Z = self.calcZ("ocean")
        TFs = {}
        
        omega = 2*C.pi*self.freqs
        sigma_s = 1/self.ocean_model["rho"]
        k2 = 1.j*omega*C.mu_0*sigma_s
        k = np.sqrt(k2)
        kd = k*self.ocean_model["depth"]
        
        for kind in kinds:
            TFs[kind] = self.functions[kind](Z, Zd, kd)
        return TFs
    
    def summary_plots(self):
        """
        Create a summary plot of the analysis
        """
        fig, ax = plotlib.create_pane()
        TFs = self.calcTF()
        plotlib.plotTF(TFs, self.freqs, self.ocean_model["rho"], self.ocean_model["depth"])
        return