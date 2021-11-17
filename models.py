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
import datetime as dt
from scipy import constants as C
import pandas as pd

from multiprocessing import Pool
from netCDF4 import Dataset
from scipy import interpolate
import glob
from decimal import Decimal

import plotlib

def toBEZpy(base="data/OceanModels/"):
    
    def fexp(number):
        (sign, digits, exponent) = Decimal(number).as_tuple()
        return len(digits) + exponent - 1
    
    def fman(number):
        return Decimal(number).scaleb(-fexp(number)).normalize()
    
    def sign(number):
        (sign, digits, exponent) = Decimal(number).as_tuple()
        return "+" if sign==0 else "-"
    
    files = glob.glob(base+"*.csv")
    files.sort()
    header = "* Lines starting with * are just comments.\n"+\
                "* Text after the numbers is ignored \n"+\
                "* BM ocean conductivity model\n"+\
                "*/ %s, INF,/ ! layer thicknesses in m\n"+\
                "*/ %s ,/ !Resistivities in Ohm-m\n"+\
                "%d                             Number of layers from surface\n"
    eachline = "\n%.7f                      Conductivity in S/m (layer %d)"+\
                "\n%.3fe%s%02d                      Layer thickness in m (layer %d)\n"
    lastline = "\n1.1220100                      Semi-infinite earth conductivity"
    for f in files:
        o = pd.read_csv(f)
        thks, rhos = ", ".join([str(x) for x in o["Thk(km)"]]), ", ".join([str(x) for x in o["Rho(ohm-m)"]])
        rhos += (", %.3f"%(1./1.1220100))
        body = header%(thks, rhos, len(o))
        for i, row in o.iterrows():
            th = row["Thk(km)"]*1e3
            body += eachline%(1/row["Rho(ohm-m)"], i+1, fman(th), sign(th), fexp(th), i+1)
            #print(row["Rho(ohm-m)"], row["Thk(km)"])
            pass
        body += lastline
        with open(f.replace(".csv", ".txt"), "w") as f: f.writelines(body)
    return

class OceanModel(object):
    """
    This class is to emulate electric and magnetic fields
    in a 1D ocean. We provide the 1D Earth conductivity model and
    one layer of ocean parameters (depth and resistivity) to emulate 
    1D ocean.
    """
    
    def __init__(self, site=None, model_name="BME", ocean_model={"depth":5e3, "rho":0.25}, flim=[1e-6, 1e0]):
        self.model_name = model_name
        self.ocean_model = ocean_model
        self.site = site if site else bezpy.mt.read_1d_usgs_profile("data/ocean_model_%s.txt"%model_name)
        self.freqs = np.linspace(flim[0], flim[1], int(flim[1]/flim[0])+1)
        self.functions = {
            "Ef2Bs": lambda Z, Zd, kd: Zd/(np.cosh(kd) + (Zd*np.sinh(kd)/Z)),
            "Bf2Bs": lambda Z, Zd, kd: 1./(np.cosh(kd) + (Zd*np.sinh(kd)/Z))
        }
        self.summary_plots()
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
    
    def get_TFs(self, key="Ef2Bs"):
        TFs = self.calcTF()
        tf = pd.DataFrame()
        tf["freq"], tf[key] = self.freqs, TFs[key]
        return tf
    
    def summary_plots(self):
        """
        Create a summary plot of the analysis
        """
        fig, ax = plotlib.create_pane()
        TFs = self.calcTF()
        plotlib.plotTF(ax, TFs, self.freqs, self.ocean_model["rho"], self.ocean_model["depth"])
        fig.savefig("docs/TFs.png", bbox_inches="tight")
        return
    
    @staticmethod
    def getOceanModel(bin_num=1):
        site = bezpy.mt.read_1d_usgs_profile("data/OceanModels/RhoZ_Bin%d.txt"%bin_num)
        ocean_model = {"depth":site.thicknesses[0], "rho":site.resistivities[0]}
        site.thicknesses, site.resistivities = site.thicknesses[1:], site.resistivities[1:]
        om = OceanModel(site=site, ocean_model=ocean_model)
        return om
    
def intp_data(x, y, xnew):
    f = interpolate.interp1d(x, y, fill_value = "extrapolate")
    return f(xnew)

from scipy.fft import fft, fftfreq, ifft
def convolve(signal_a_time, signal_b_freq, b_freq):
    signal_a_freq = fft(signal_a_time)
    frq = fftfreq(len(signal_a_time), 1)
    for bf, f in zip(signal_b_freq, b_freq):
        i = np.argmin(abs(f-frq))
        signal_a_freq[i] = signal_a_freq[i]*bf
    signal_ab_time = ifft(signal_a_freq)
    return signal_ab_time

class BFieldAnalysis(object):
    """
    Analyze B-fields for various radar stations
    1. Read data using bezpy.read_iaga
    2. Convert WDC to XYZ format
    3. Plot as a stack
    4. Create FFT
    """
    
    def __init__(self, stns=["FRD", "STJ", "HAD"]):
        """
        Initialize the parameters
        """
        self.stns = stns
        
        self.files = dict(
            FRD = ["FRD_19890312_XYZ.txt", "FRD_19890313_XYZ.txt", "FRD_19890314_XYZ.txt"],
            STJ = ["STJ_19890312_XYZ.txt", "STJ_19890313_XYZ.txt", "STJ_19890314_XYZ.txt"],
            HAD = ["HAD_19890312_HDZ.txt", "HAD_19890313_HDZ.txt", "HAD_19890314_HDZ.txt"],
        )
        self.frames = dict()
        self.read_dataset()
        #self.fft()
        return
    
    def read_dataset(self):
        """
        Read dataset 
        """
        base = "data/1989/"
        for stn in self.stns:
            o = pd.DataFrame()
            for f in self.files[stn]:
                fname = base + f
                if os.path.exists(fname): o = pd.concat([o, bezpy.mag.read_iaga(fname)])
            self.frames[stn] = o
        self.stack_plots()
        return
    
    def stack_plots(self, kind="Bxy"):
        if kind == "Bxy": plotlib.plot_xy_magnetic_field_oneplot(self.stns, self.frames)
        if kind == "Bxyfft": plotlib.plot_xy_magnetic_field_fft(self.stns, self.fft_frames)
        if kind == "Exyfft": plotlib.plot_xy_electric_field_fft(self.stns, self.Efft_frames)
        if kind == "Exy": plotlib.plot_xy_electric_field(self.stns, self.E_frames)
        return
    
    def fft(self):
        self.fft_frames = {}
        for stn in self.stns:
            self.fft_frames[stn] = {}
            for a in ["X", "Y"]:
                sig = np.array(self.frames[stn][a])
                N = len(sig)
                Xf = np.fft.fftshift(np.fft.fft(sig))
                Fs = 1/60. # in Hz
                freq = np.arange(-Fs/2, Fs/2, Fs/N)
                self.fft_frames[stn][a] = {"Fs": Fs, "Xf": Xf[N//2:], "freq": freq[N//2:], "N": N}
        self.stack_plots(kind="Bxyfft")
        return
    
    def calculate_Ef(self, tf, key="Ef2Bs", flim=[1e-6, 8e-3]):
        freqs = np.linspace(flim[0], flim[1], int(flim[1]/flim[0])+1)
        tf = tf[(tf.freq>=flim[0]) & (tf.freq<=flim[1])]
        ntf = pd.DataFrame()
        ntf["freq"], ntf[key] = freqs, intp_data(np.array(tf.freq), np.array(tf[key]), freqs)
        self.E_frames = {}
        for stn in self.stns:
            frame, self.E_frames[stn] = self.frames[stn], {}
            ufrm = pd.DataFrame()
            ufrm["Time"] = frame.index.tolist()
            for a in ["X", "Y"]:
                ufrm[a] = convolve(np.array(frame[a]), ntf[key], ntf.freq)
            ufrm = ufrm.set_index("Time")
            self.E_frames[stn] = ufrm
        self.stack_plots(kind="Exy")
        return
    
if __name__ == "__main__":
    #om = OceanModel()
    #BFieldAnalysis().calculate_Ef(om.get_TFs())
    OceanModel.getOceanModel(2)
    pass