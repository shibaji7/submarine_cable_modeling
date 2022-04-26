"""
    synthetic.py: Module is used to implement Synthetic B- and E-Field structures.
    Class:
    ------
    SynB: Synthetic B-Field generations
    SynE: Synthetic E-Field generations
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
import sys
sys.path.extend(["py/", "py/config/"])
# Import required packages
import bezpy
import numpy as np
from scipy import constants as C
import datetime as dt
import pandas as pd
from scipy.stats import pearsonr
from loguru import logger

import utility
from oml import OceanModel
from plotlib import BfieldSummary



class SynB(object):
    """
    This class is dedicated to synthetic B-field simulations and 
    storing the results.
    """
    
    def __init__(self, args):
        self.args = args
        logger.info(f"Synthetic B-field run parameters")
        return
    
    def create_synthetic_B_field(self, Am, Phim, Tm, t):
        """
        This function is responsible for creating a 
        synthetic magnetic (B) field with following sets 
        of parameters.

        Parameter:
        ----------
        Am (list) - Magntitue at different freuency components (m)
        Phim (list) - Phase at different freuency components (m)
        Tm (list) - Periods of different freuency components (m)
        t (float) - Total time
        """
        self.Bfield = {}
        self.stations = ["syn"]
        self.Am, self.Phim, self.Tm = Am, Phim, Tm
        t = np.linspace(0, t, t,endpoint=False)
        Bt = np.zeros(len(t))
        for A, Phi, T in zip(Am, Phim, Tm):
            Bt += A*np.sin(2*np.pi*t/T + np.deg2rad(Phi))
        self.Bfield["syn"] = pd.DataFrame()
        self.Bfield["syn"]["X"], self.Bfield["syn"]["dTime"] = Bt, t
        self.components = ["X"]
        self.w = utility.get_tapering_function(np.array(self.Bfield["syn"]["dTime"]), 
                                               self.args.Bfield.tapering)
        return
    
    def read_dataset(self):
        """
        Read all the files from the dataset files
        """
        self.Bfield = {}
        sources = self.args.Bfield.data_sources.sources
        self.stations = self.args.Bfield.data_sources.stns
        for f in sources:
            stn = f.split("/")[-1].split("_")[0]
            if stn not in self.Bfield.keys(): self.Bfield[stn] = pd.DataFrame()
            if os.path.exists(f): self.Bfield[stn] = pd.concat([self.Bfield[stn], bezpy.mag.read_iaga(f)])
        self.components = ["X", "Y"]
        for stn in self.stations:
            Bf = self.Bfield[stn]
            t0 = Bf.index.tolist()[0].to_pydatetime()
            Bf["dTime"] = [(x.to_pydatetime()-t0).total_seconds() 
                           for x in Bf.reset_index().Time.tolist()]
            self.w = utility.get_tapering_function(np.array(Bf["dTime"]), 
                                                   self.args.Bfield.tapering)
            T0 = np.array(Bf["dTime"])
            for c in self.components:
                Bf[c] = utility.detrend_magnetic_field(np.array(Bf[c]), T0,
                                                       self.args.Bfield.tapering)
            self.Bfield[stn] = Bf
        return
    
    def run(self):
        """
        Run all the steps
        """
        self.syn = False
        # Genarting / reading all B-field dataset
        if hasattr(self.args.Bfield, "structure"):
            self.syn = True
            Am, Tm, Phim, t = np.array(self.args.Bfield.structure.Am),\
                    np.array(self.args.Bfield.structure.Tm_min)*60,\
                    np.array(self.args.Bfield.structure.Phim),\
                    int(np.array(self.args.Bfield.structure.T_hours)*60*60)
            self.create_synthetic_B_field(Am, Phim, Tm, t)
        else: self.read_dataset()
        bdir = self.args.out_dir
        # Process B-field data
        if hasattr(self.args.Bfield, "models"):
            base_model_path = self.args.Bfield.models.model_dir_path
            self.Efield = {}
            for model in self.args.Bfield.models.model_list: # run for each model
                d, r, m = model.ocean_depth, model.ocean_resistivity, model.earth_model
                prep = "h-%d.r-%.2f"%(d, r)
                logger.info(f"Model params:{prep}")
                mpath = base_model_path%"BME" if "uniform" in m else base_model_path%m
                om = OceanModel(model.earth_model, mpath, ocean_model={"depth":d, "rho":r})
                if "uniform" in m:
                    rho = float(m.split(".")[-1])
                    om.site.resistivities[:] = rho
                logger.info(f"Synthetic B {m}->OM({om.model_name})")
                for stn in self.stations:
                    Enum = self.solve_numerical_Et(om, stn)
                    if self.syn:
                        km, tf = self.draw_Km_table(om), self.draw_TF_table(om)
                        Eanl = self.solve_analytical_Et(tf, stn)
                        r = self.check_analytical_numerical(Eanl, Enum)
                        pname = bdir + "summary_plot_%s_%s_%s.png"%(m, prep, stn)
                        self.summary_plot(om, Eanl, Enum, pname, stn)
                        km.to_csv(bdir + "Kfm_%s_%s_%s.csv"%(m, prep, stn), float_format="%g", 
                                  header=True, index=False)
                        tf.to_csv(bdir + "TF_%s_%s_%s.csv"%(m, prep, stn), float_format="%g", 
                                  header=True, index=False)
                        Eanl.to_csv(bdir + "Eanl_%s_%s_%s.csv"%(m, prep, stn), float_format="%g",
                                    header=True, index=False)
                        np.savetxt(bdir + "cor_%s_%s_%s.csv"%(m, prep, stn), np.array([r]).T, 
                               delimiter=",", header="rho", comments="")
                    Enum.to_csv(bdir + "Enum_%s_%s_%s.csv"%(m, prep, stn), float_format="%g", 
                                header=True, index=False)
                    self.Bfield[stn].to_csv(bdir + "Bt_%s.csv"%stn, float_format="%g", header=True, index=False)
                    np.savetxt(bdir + "Wp_%s.csv"%stn, np.array([np.array(self.Bfield[stn].dTime), self.w]).T, 
                               fmt="%g", delimiter=",", header="t,Wp", comments="")
        else:
            for stn in self.stations:
                self.Bfield[stn].to_csv(bdir + "Bt_%s.csv"%stn, float_format="%g", header=True, index=False)
                np.savetxt(bdir + "Wp_%s.csv"%stn, np.array([np.array(self.Bfield[stn].dTime), self.w]).T, 
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
        o["Amplitude (mV/km/nT)"], o["Phase (deg)"] = np.absolute(o.E2B),\
                        np.angle(o.E2B, deg=True)
        o = o.rename(columns={"freq":"fm (Hz)"}).drop(columns=["E2B"])
        o.index.name = "m"
        return o.copy()
    
    def solve_analytical_Et(self, tf, stn):
        """
        Solve for analytical Et from Bt and TF (tf)
        """
        Et = np.zeros(len(self.Bfield[stn]["X"]))
        m = 0
        t = np.array(self.Bfield[stn].dTime)
        for A, Phi, T in zip(self.Am, self.Phim, self.Tm):
            Et += np.absolute(tf["Amplitude (mV/km/nT)"])[m]*A*\
                    np.sin(2*np.pi*t/T + np.deg2rad(Phi + np.deg2rad(tf["Phase (deg)"])[m]))
            m += 1
        o = pd.DataFrame()
        o["dTime"], o["X"] = self.Bfield[stn].dTime, Et
        return o.copy()
    
    def solve_numerical_Et(self, om, stn):
        """
        Solve for Et from FFT and IFFT operations
        """
        o = pd.DataFrame()
        o["dTime"] = self.Bfield[stn].dTime
        for c in self.components:
            Bt, dT = np.array(self.Bfield[stn][c]) * self.w,\
                self.Bfield[stn].dTime.iloc[1]-self.Bfield[stn].dTime.iloc[0]
            Bf, f = utility.fft(Bt, dT)
            E2B = np.array(om.get_TFs(freqs=f).E2B)
            Et = utility.ifft(E2B*Bf)
            o[c] = Et
        return o.copy()
    
    def check_analytical_numerical(self, Eanl, Enum):
        """
        Regression check for numerical and 
        """
        L = int(len(Eanl)/3)
        for c in self.components:
            r, _ = pearsonr(Eanl[c].tolist()[L:-L], Enum[c].tolist()[L:-L])
            logger.info(f"Corr(Eanl,Enum){'%.10f'%r})")
        return r
    
    def summary_plot(self, om, Eanl, Enum, pname, stn):
        """
        Create summary plots
        """
        summ = BfieldSummary()
        ax = summ.add_Bfield_Seq(self.Bfield[stn], Enum)
        ax.text(0.01, 1.05, "Numerical Sol.", ha="left", va="center", transform=ax.transAxes)
        ax = summ.add_Bfield_Seq(self.Bfield[stn], Eanl)
        ax.text(0.01, 1.05, "Analytical Sol.", ha="left", va="center", transform=ax.transAxes)
        summ.add_Es(Eanl.X.tolist()[1000:-1000], Eanl.X.tolist()[1000:-1000])
        summ.add_TK_param(om.get_TFs())
        summ.save(pname)
        summ.close()
        return
    
class SynE(object):
    """
    This class is dedicated to synthetic E-field simulations and 
    storing the results.
    """
    
    def __init__(self, args, Efield=None):
        self.args = args
        self.Efield = Efield
        logger.info(f"Synthetic E-field run parameters")
        return
    
    def create_synthetic_E_field(self):
        """
        create synthetic E field
        """
        Eframe = pd.DataFrame()
        Eframe["X"], Eframe["Y"], Eframe["Time"] = self.args.Efield.structure.X,\
                self.args.Efield.structure.Y, [dt.datetime.now()]*len(self.args.Efield.structure.X)
        Eframe = Eframe.set_index("Time")
        Eframe["dTime"] = np.arange(len(self.args.Efield.structure.X))
        self.Efield = { "syn": Eframe }
        self.components = ["X", "Y"]
        self.stations = [ "syn" ]
        return
    
    def read_dataset(self):
        """
        Read all the files from the dataset files
        """
        self.Efield = {}
        sources = self.args.Bfield.data_sources.sources
        self.stations = self.args.Bfield.data_sources.stns
        self.components = ["X", "Y"]
        for f in sources:
            stn = f.split("/")[-1].split("_")[0]
            if stn not in self.Efield.keys(): self.Efield[stn] = pd.DataFrame()
            if os.path.exists(f): 
                self.Efield[stn] = pd.concat([self.Efield, 
                                              pd.read_csv(f, index=["Time"], parse_dates=["Time"])])
        for stn in self.stations:
            t0 = self.Efield[stn].index.tolist()[0].to_pydatetime()
            self.Efield[stn]["dTime"] = [(x.to_pydatetime()-t0).total_seconds() 
                                         for x in self.Efield[stn].reset_index().Time.tolist()]
        return
    
    def run(self):
        """
        Run all the steps
        """
        # Genarting / reading all E-field dataset
        if hasattr(self.args.Efield, "structure"):
            self.syn = True
            self.create_synthetic_E_field()
        else: self.read_dataset()
        bdir = self.args.out_dir
        for stn in self.stations:
            self.Efield[stn].to_csv(bdir + "Efield_%s.csv"%stn, float_format="%g", header=True, index=False)
        return