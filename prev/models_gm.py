"""models.py: Module is used to implement various sea-earth models and calculate TFs"""

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

model_map = dict(
        CS_W = 1,
        DO_1 = 4,
        DO_2 = 3,
        DO_3 = 2,
        DO_4 = 5,
        MAR = 7,
        CS_E = 9
        )

bin_to_station_map = {1: "FRD", 2: "FRD", 3: "STJ",
                      4: "STJ", 5: "STJ", 6: "STJ",
                      7: "STJ", 8: "HAD", 9: "HAD"}
bin_eastern_edge = {1: (38.79, -72.62), 2: (37.11, -68.94), 3: (39.80, -48.20), 
                    4: (40.81, -45.19), 5: (43.15, -39.16), 6: (44.83, -34.48), 
                    7: (46.51, -22.43), 8: (47.85, -9.05), 9: (50.79, -4.55)}
bin_western_edge = {1: (39.6, -74.33), 2: (38.79, -72.62), 3: (37.11, -68.94), 
                    4: (39.80, -48.20), 5: (40.81, -45.19), 6: (43.15, -39.16), 
                    7: (44.83, -34.48), 8: (46.51, -22.43), 9: (47.85, -9.05)}
bin_ocean_depth = {1: 100, 2: 4000, 3: 5200,
                   4: 4000, 5: 4800, 6: 4000, 
                   7: 3000, 8: 4500, 9: 100}

def segmentL(binsec=1):
    """
    Northward (Ln) and Eastward (Le) distance (km) between the
    ends of the cable section
    """
    # PARAMETERS OF THE WGS84 EARTH MODEL
    a = 6378.137 # Equatorial radius
    b = 6356.752 # Polar Radius
    e = np.sqrt(0.00669437999014) # Eccentricity
    east_edge, west_edge = bin_eastern_edge[binsec], bin_western_edge[binsec]
    phi = 0.5*(west_edge[1]+east_edge[1])
    Ln = (111.133-0.56*np.cos(np.deg2rad(2*phi)))*np.abs(east_edge[0]-west_edge[0])
    Le = (111.5065-0.1872*np.cos(np.deg2rad(2*phi)))*np.cos(np.deg2rad(phi))*np.abs(east_edge[1]-west_edge[1])
    return Ln, Le

def toBEZpy(base="data/OceanModels/"):
    
    def fexp(number):
        (sign, digits, exponent) = Decimal(number).as_tuple()
        return len(digits) + exponent - 1
    
    def fman(number):
        return Decimal(number).scaleb(-fexp(number)).normalize()
    
    def sign(number):
        (sign, digits, exponent) = Decimal(number).as_tuple()
        return "+" if sign==0 else "-"
    
    files = glob.glob(base+"*Bin*.csv")
    files.sort()
    header = "* Lines starting with * are just comments.\n"+\
                "* Text after the numbers is ignored \n"+\
                "* BM ocean conductivity model\n"+\
                "*/ %s, INF,/ ! layer thicknesses in km\n"+\
                "*/ %s ,/ !Resistivities in Ohm-m\n"+\
                "%d                             Number of layers from surface\n"
    eachline = "\n%.7f                      Conductivity in S/m (layer %d)"+\
                "\n%.3fe%s%02d                      Layer thickness in m (layer %d)\n"
    lastline = "\n1.1220100                      Semi-infinite earth conductivity"
    ocean_layer = []
    for f in files:
        o = pd.read_csv(f)
        bname = f.split("_")[-1].replace(".csv","")
        ocean_layer.append({"bin": bname, "depth": o["Thk(km)"][0]*1e3, "rho": o["Rho(ohm-m)"][0]})
        thks, rhos = ", ".join([str(x) for x in o["Thk(km)"][1:]]),\
                ", ".join([str(x) for x in o["Rho(ohm-m)"][1:]])
        rhos += (", %.3f"%(1./1.1220100))
        body = header%(thks, rhos, (len(o)-1))
        for i, row in o.iterrows():
            if i > 0:
                th = row["Thk(km)"]*1e3
                body += eachline%(1/row["Rho(ohm-m)"], i, fman(th), sign(th), fexp(th), i)
        body += lastline
        with open(f.replace(".csv", ".txt"), "w") as f: f.writelines(body)
    ocean_layer = pd.DataFrame.from_records(ocean_layer)
    ocean_layer.to_csv("data/OceanModels/OceanLayers.csv", header=True, index=False)
    return

class OceanModel(object):
    """
    This class is to emulate electric and magnetic fields
    in a 1D ocean. We provide the 1D Earth conductivity model and
    one layer of ocean parameters (depth and resistivity) to emulate 
    1D ocean.
    """
    
    def __init__(self, site=None, model_name="BME", 
                 ocean_model={"depth":5e3, "rho":0.25}, 
                 flim=[1e-4, 1e-2]):
        self.model_name = model_name
        self.ocean_model = ocean_model
        self.site = bezpy.mt.read_1d_usgs_profile("data/ocean_model_%s.txt"%model_name) if site is None else site
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
        ocean_info = pd.read_csv("data/OceanModels/OceanLayers.csv")
        ocean_info = ocean_info[ocean_info.bin=="Bin"+str(bin_num)]
        site = bezpy.mt.read_1d_usgs_profile("data/OceanModels/RhoZ_Bin%d.txt"%bin_num)
        ocean_model = {"depth":ocean_info.depth.tolist()[0], "rho":ocean_info.rho.tolist()[0]}
        om = OceanModel(site=site, ocean_model=ocean_model)
        return om


class BFieldAnalysis(object):
    """
    Analyze B-fields for various radar stations
    1. Read data using bezpy.read_iaga
    2. Convert WDC to XYZ format
    3. Plot as a stack
    4. Compute Et
    """
    
    def __init__(self, dates, bin_id, db="data/{year}/", coord="XYZ", plot=True):
        """
        Initialize the parameters
        """
        self.bin_id = bin_id
        self.stn = bin_to_station_map[self.bin_id]
        self.files = []
        self.plot = plot
        for d in dates:
            f = db.format(year=d.year) + "{stn}_{date}_{coord}.txt".format(stn=self.stn, 
                                                                           date=d.strftime("%Y%m%d"), 
                                                                           coord=coord)
            self.files.append(f)   
        self.B_frame = None
        self.read_dataset()
        return
    
    def read_dataset(self):
        """
        Read dataset 
        """
        o = pd.DataFrame()
        for f in self.files:
            if os.path.exists(f): o = pd.concat([o, bezpy.mag.read_iaga(f)])
        self.B_frame = o
        if self.plot: self.stack_plots()
        self.detrain_B()
        return
    
    def stack_plots(self, kind="Bxy"):
        if kind == "Bxy": plotlib.plot_Bxy(self.stn, self.B_frame, fname="docs/Bxy.Field.B%02d.png"%self.bin_id)
        if kind == "Exy": plotlib.plot_Exy(self.stn, self.E_frame, fname="docs/Exy.Field.B%02d.png"%self.bin_id)
        if kind == "BExy": plotlib.plot_BExy(self.stn, self.B_frame, self.E_frame, 
                                            fname="docs/BExy.Field.B%02d.png"%self.bin_id)
        return
    
    def compute_Et(self):
        """
        Compute Et via RFFT and IRFFT block 
        om: Ocean Model
        """
        self.om = OceanModel.getOceanModel(self.bin_id)
        o = pd.DataFrame()
        n = len(self.B_frame)
        dT = (self.B_frame.index.tolist()[1]-self.B_frame.index.tolist()[0]).total_seconds()
        for a in ["X", "Y"]:
            Baf = 2.0/n * np.fft.rfft(self.B_frame[a])
            frq = np.fft.rfftfreq(len(self.B_frame[a]))/dT
            frq[0] = frq[1]
            #print(f" Frequency range of the fields: {np.min(frq)}~{np.max(frq)}; dF:{1./dT}")
            tx = self.om.get_TFs(freqs=frq)
            Eaf = Baf*np.array(tx.Ef2Bs)
            #print(f" TF abs range: {np.min(np.abs(tx.Ef2Bs))}~{np.max(np.abs(tx.Ef2Bs))}")
            #print(f" Bf abs range: {np.min(np.abs(Baf))}~{np.max(np.abs(Baf))}")
            #print(f" Ef abs range: {np.min(np.abs(Eaf))}~{np.max(np.abs(Eaf))}")
            o[a] = np.fft.irfft(Eaf)*n/2
        o["Time"] = self.B_frame.index.tolist()
        self.E_frame = o.set_index("Time")
        if self.plot: self.stack_plots("BExy")
        self.compute_Vj()
        return
    
    def compute_Vj(self):
        """
        Calculate total electric potential induced along the cable segment.
        Vj = Ej_n(t)Lj_n + Ej_e(t)Lj_e
        """
        Ln, Le = segmentL(self.bin_id)
        self.Vj = np.array(self.E_frame.X)*Ln + np.array(self.E_frame.Y)*Le
        self.E_frame["Vj"] = self.Vj
        return
    
    def detrain_B(self, p=0.1):
        """
        Preprocess the B-field before the analysis
        1. Remove mean
        2. Detrain
        3. Tapering the ends for removing sporious frequency leakage
        """
        start = self.B_frame.index.tolist()[0]
        self.B_frame = self.B_frame.reset_index()
        self.B_frame["Tm"] = self.B_frame.Time.apply(lambda x: (x-start).total_seconds())
        tm = np.array(self.B_frame["Tm"])
        self.B_frame = self.B_frame.set_index("Time")
        T = len(self.B_frame)
        P, P2 = int(T*p), int(T*p/2)
        wp = np.zeros(T)
        wp[:P2] = 0.5*(1 - np.cos(2*np.pi*tm[:P2]/P))
        wp[P2:T-P2] = 1.
        wp[T-P2:] = 0.5*(1 - np.cos(2*np.pi*(tm[-1]-tm[T-P2:])/P))
        for a in ["X", "Y"]:
            dat = self.B_frame[a]
            dat = dat - np.median(dat[:120])
            self.B_frame[a] = dat*wp
        return
    
if __name__ == "__main__":
    #toBEZpy()
    for b, c in zip(range(1,10), ["XYZ"]*7 + ["HDZ"]*2):
        BFieldAnalysis([dt.datetime(1989,3,12),dt.datetime(1989,3,13),dt.datetime(1989,3,14)],
                       b, coord=c).compute_Et()
        break