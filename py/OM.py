"""OM.py: Module is used to implement Ocean Model with Layered-Earth structures"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import sys
sys.path.extend(["py/", "py/config/"])

# Import required packages
import bezpy
import numpy as np
from scipy import constants as C
import pandas as pd

from scipy.stats import pearsonr
from loguru import logger

import utility

class OceanModel(object):
    """
    This class is to emulate electric and magnetic fields
    in a 1D ocean. We provide the 1D Earth conductivity model and
    one layer of ocean parameters (depth and resistivity) to emulate 
    1D ocean.
    """
    
    def __init__(self, model_name, model_path, ocean_model={"depth":5e3, "rho":0.25}, 
                 site=None, flim=[1e-4, 1e-2]):
        self.model_name = model_name
        self.ocean_model = ocean_model
        self.site = bezpy.mt.read_1d_usgs_profile(model_path) if site is None else site
        self.freqs = np.linspace(flim[0], flim[1], int(flim[1]/flim[0])+1)
        self.functions = {
            "Ef2Bs": lambda Z, Zd, kd: Zd/(np.cosh(kd) + (Zd*np.sinh(kd)/Z)),
            "Bf2Bs": lambda Z, Zd, kd: 1./(np.cosh(kd) + (Zd*np.sinh(kd)/Z))
        }
        return
    
    def calcZ(self, layer="ocean", freqs=None):
        freqs = np.copy(self.freqs) if freqs is None else freqs
        omega = 2*C.pi*freqs
        sigma_s = 1/self.ocean_model["rho"]
        k2 = 1.j*omega*C.mu_0*sigma_s
        k = np.sqrt(k2)
        Zo = (1.j*omega*C.mu_0/k)/(C.mu_0/1.e-3)
        Kf = self.site.calcZ(freqs)[1]
        self.Z = {"ocean": Zo, "floor": Kf}
        return self.Z[layer]
    
    def calcTF(self, kinds=["Ef2Bs", "Bf2Bs"], freqs=None):
        freqs = np.copy(self.freqs) if freqs is None else freqs
        omega = 2*C.pi*freqs
        Zd = self.calcZ("floor", freqs=freqs)
        Z = self.calcZ("ocean", freqs=freqs)
        TFs = {}
        
        omega = 2*C.pi*freqs
        sigma_s = 1/self.ocean_model["rho"]
        k2 = 1.j*omega*C.mu_0*sigma_s
        k = np.sqrt(k2)
        kd = k*self.ocean_model["depth"]
        
        for kind in kinds:
            TFs[kind] = self.functions[kind](Z, Zd, kd)
        return TFs
    
    def get_TFs(self, key="Ef2Bs", freqs=None):
        TFs = self.calcTF(freqs=freqs)
        tf = pd.DataFrame()
        tf["freq"], tf[key] = np.copy(self.freqs) if freqs is None else freqs, TFs[key]
        return tf
    

class SynB(object):
    """
    This class is dedicated to synthetic B-field simulations and 
    storing the results.
    """
    
    def __init__(self, o, tp, v, bdir):
        self.tp = tp
        self.v = v
        self.bdir = bdir
        if self.v: logger.info(f"Synthetic B-field run parameters")
        self = utility.set_dict(self, o)
        if self.v:
            for k in vars(self).keys():
                print("     ", k, "->", vars(self)[k])
        return
    
    def run(self):
        """
        Run all the steps
        """
        self.Am, self.Tm, self.Phim = np.array(self.Am), np.array(self.Tm_min)*60, np.array(self.Phim)
        self.Bt, self.t = utility.create_synthetic_B_field(self.Am, self.Phim, self.Tm)
        self.w = utility.get_tapering_function(self.t)
        for m in self.earth_models: # run for each earth model
            for d in self.ocean_depths: # run for each ocean depth 
                for r in self.ocean_resistivities: # run for each ocean resistivity
                    prep = "h-%d.r-%.2f"%(d,r)
                    mpath = self.model_path%"BME" if "uniform" in m else self.model_path%m
                    om = OceanModel(m, mpath, ocean_model={"depth":d, "rho":r})
                    if "uniform" in m:
                        rho = float(m.split(".")[-1])
                        om.site.resistivities[:] = rho
                    if self.v: logger.info(f"Synthetic B {m}->OM({om.model_name})")
                    km, tf = self.draw_Km_table(om), self.draw_TF_table(om)
                    Eanl = self.solve_analytical_Et(tf)
                    Enum = self.solve_numerical_Et(om)
                    r = self.check_analytical_numerical(Eanl, Enum)
                    if self.save: 
                        km.to_csv(self.bdir + "Kfm_%s_%s.csv"%(m, prep), float_format="%g")
                        tf.to_csv(self.bdir + "TF_%s_%s.csv"%(m, prep), float_format="%g")
                        Eanl.to_csv(self.bdir + "Eanl_%s_%s.csv"%(m, prep), float_format="%g")
                        Enum.to_csv(self.bdir + "Enum_%s_%s.csv"%(m, prep), float_format="%g")
                        np.savetxt(self.bdir + "cor_%s_%s.csv"%(m, prep), np.array([r]).T, 
                                   delimiter=",", header="rho", comments="")
        if self.save:
            np.savetxt(self.bdir + "Bt.csv", np.array([self.t, self.Bt]).T, 
                       fmt="%g", delimiter=",", header="t,Bt", comments="")
            np.savetxt(self.bdir + "Wp.csv", np.array([self.t, self.w]).T, 
                       fmt="%g", delimiter=",", header="t,Wp", comments="")
        return
    
    def draw_Km_table(self, qx):
        """
        Create the Km(f) tables at specific frequencies.
        """
        fm = 1./self.Tm
        Kfm = qx.calcZ(layer="floor", freqs=fm)
        o = pd.DataFrame()
        o["fm (Hz)"], o["Amplitude |Km| (mV/km/nT)"], o["Phase (deg)"] = fm, np.absolute(Kfm),\
                            np.angle(Kfm, deg=True)
        o.index.name = "m"
        return o.copy()
    
    def draw_TF_table(self, qx):
        """
        Create the TF(f) tables at specific frequencies.
        """
        fm = 1./self.Tm
        o = qx.get_TFs(freqs=fm)
        o["Amplitude (mV/km/nT)"], o["Phase (deg)"] = np.absolute(o.Ef2Bs),\
                        np.angle(o.Ef2Bs, deg=True)
        o = o.rename(columns={"freq":"fm (Hz)"}).drop(columns=["Ef2Bs"])
        o.index.name = "m"
        return o.copy()
    
    def solve_analytical_Et(self, tf):
        """
        Solve for analytical Et from Bt and TF (tf)
        """
        Et = np.zeros(len(self.t))
        m = 0
        for A, Phi, T in zip(self.Am, self.Phim, self.Tm):
            Et += np.absolute(tf["Amplitude (mV/km/nT)"])[m]*A*\
                    np.sin(2*np.pi*self.t/T + np.deg2rad(Phi + np.deg2rad(tf["Phase (deg)"])[m]))
            m += 1
        o = pd.DataFrame()
        o["t"], o["Et"] = self.t, Et
        return o.copy()
    
    def solve_numerical_Et(self, om):
        """
        Solve for Et from FFT and IFFT operations
        """
        Bt, dT = self.Bt * self.w, self.t[1]-self.t[0]
        Bf, f = utility.fft(Bt, dT)
        Ef2Bs = np.array(om.get_TFs(freqs=f).Ef2Bs)
        Et = utility.ifft(Ef2Bs*Bf)
        o = pd.DataFrame()
        o["t"], o["Et"] = self.t, Et
        return o.copy()
    
    def check_analytical_numerical(self, Eanl, Enum):
        """
        Regression check for numerical and 
        """
        L = int(len(Eanl)/3)
        r, _ = pearsonr(Eanl.Et.tolist()[L:-L], Enum.Et.tolist()[L:-L])
        if self.v: logger.info(f"Corr(Eanl,Enum){'%.10f'%r})")
        return r