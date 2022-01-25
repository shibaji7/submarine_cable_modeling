"""cable.py: Module is used to implement cable section analysis and an event study"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
import pandas as pd
import numpy as np
from loguru import logger

import bezpy

from oml import OceanModel
from dstl import TransmissionLine, NodalAnalysis
import utility
import plotlib

class CableSection(object):
    """
    This Class is dedicated to cable section analysis including
    1. B-field data load.
    2. Electric field calculation
    3. Induced voltage calculation
    4. DSTL voltage calculation
    5. Total voltage drop calculation
    """
    
    def __init__(self, model_num, model_path, stn_name, ocean_model,
                 east_edge, west_edge, data_files):
        self.model_num = model_num
        self.model_path = model_path
        self.stn_name = stn_name
        self.ocean_model = ocean_model
        self.east_edge = east_edge
        self.west_edge = west_edge
        self.data_files = data_files
        self.load_om()
        self.load_dataset()
        self.components = ["X", "Y"]
        return
    
    def load_om(self, kind="floor"):
        """
        Load Ocean Model with base parameters
        """
        name = "Bin%02d"%self.model_num
        self.om = OceanModel(name, self.model_path, ocean_model=self.ocean_model, kind=kind)
        return
    
    def load_dataset(self):
        """
        Load dataset from given data files
        """
        self.B = pd.DataFrame()
        for f in self.data_files:
            if os.path.exists(f): self.B = pd.concat([self.B, bezpy.mag.read_iaga(f)])
        return
    
    def compute_Vj(self):
        """
        Calculate total electric potential induced along the cable segment.
        Vj = Ej_n(t)Lj_n + Ej_e(t)Lj_e
        """
        Ln, Le = self.segmentL()
        self.V = pd.DataFrame()
        self.V["Time"] = self.B.index.tolist()
        self.V["Vj"] = np.array(self.E.X)*Ln + np.array(self.E.Y)*Le
        self.V = self.V.set_index("Time")
        return
    
    def compute_numerical_Et(self):
        """
        Compute Et using numerical FFT and IFFT block 
        """
        self.E = pd.DataFrame()
        self.E["Time"] = self.B.index.tolist()
        stime = self.E.Time.tolist()[0]
        t = np.array(self.E.Time.apply(lambda x: (x-stime).total_seconds()))
        for a in self.components:
            Bt = utility.detrend_magnetic_field(np.array(self.B[a]), t)
            dT = (self.B.index.tolist()[1]-self.B.index.tolist()[0]).total_seconds()
            Bf, f = utility.fft(Bt, dT)
            E2B = np.array(self.om.get_TFs(freqs=f).E2B)
            Et = 2*utility.ifft(E2B*Bf)
            self.E[a] = Et
        self.E = self.E.set_index("Time")
        self.compute_Vj()
        return
    
    def segmentL(self):
        """
        Northward (Ln) and Eastward (Le) distance (km) between the
        ends of the cable section
        """
        # PARAMETERS OF THE WGS84 EARTH MODEL
        a = 6378.137 # Equatorial radius
        b = 6356.752 # Polar Radius
        e = np.sqrt(0.00669437999014) # Eccentricity
        phi = 0.5*(self.west_edge[1] + self.east_edge[1])
        Ln = (111.133-0.56*np.cos(np.deg2rad(2*phi)))*np.abs(self.east_edge[0]-self.west_edge[0])
        Le = (111.5065-0.1872*np.cos(np.deg2rad(2*phi)))*np.cos(np.deg2rad(phi))*\
                np.abs(self.east_edge[1]-self.west_edge[1])
        return Ln, Le

class EventAnalysis(object):
    """
    This model is dedicated to run an event analsys (e.g., 1989 Dec 13)
    """
    
    def __init__(self, o, bdir, v):
        self.bdir = bdir
        self.v = v
        self = utility.set_dict(self, o)
        if self.v: 
            logger.info(f"Event analysis and cable parameters")
            utility.print_rec(self)
        ol = pd.read_csv(self.ocean_layers_csv)
        self.tx_lines = []
        for b, depth, rho in zip(self.bins, ol.depth, ol.rho):
            om = {"depth":depth, "rho":rho}
            bd = vars(self.bin_details)[str(b)]
            model_path = self.model_location%b
            data_files = []
            for date in self.data.dates:
                date = date.split("-")
                fn = self.data.file.format(stn=bd.stn, coord=bd.coord, 
                                           y=date[0], m=date[1], d=date[2])
                data_files.append(fn)
            cs = CableSection(b, model_path, bd.stn, om, bd.eastern_edge, 
                              bd.western_edge, data_files)
            cs.compute_numerical_Et()
            setattr(vars(self.bin_details)[str(b)], "cable_sec", cs)
            tx = TransmissionLine(b, cs.om.site, bd, om)
            tx.compiles(cs.E)
            setattr(vars(self.bin_details)[str(b)], "tx_line", tx)
            self.tx_lines.append(tx)
        return
    
    def calclulate_total_parameters(self):
        """
        Method is responsible for total E-field, V, and 
        V[e-e] caclulations and nodal analysis
        """
        # Nodal analysis
        self.nodal_analysis = NodalAnalysis(self.tx_lines, self.bdir)
        self.nodal_analysis.equivalent_nodel_analysis()
        self.nodal_analysis.solve_admitance_matrix()
        self.nodal_analysis.consolidate_final_result()
        U0, U1 = self.nodal_analysis.get_voltage_ends_of_cable("Y")
        # Toral parameter calculations
        self.tot_params = pd.DataFrame()
        self.tot_params["Time"] = vars(self.bin_details)["1"].cable_sec.B.index.tolist()
        Ex, Ey, V = np.zeros(len(self.tot_params)), np.zeros(len(self.tot_params)),\
                    np.zeros(len(self.tot_params))
        for b in self.bins:
            cs = vars(self.bin_details)[str(b)].cable_sec
            Ex += np.array(cs.E.X)
            Ey += np.array(cs.E.Y)
            V += np.array(cs.V.Vj)
        self.tot_params["Ex(mv/km)"], self.tot_params["Ey(mv/km)"],\
                self.tot_params["Vc(V)"] = Ex, Ey, V/1000.
        self.tot_params["Vt(V)"] = U0-U1+(V/1000.)
        self.tot_params = self.tot_params.set_index("Time")
        if self.save: self.save_data()
        return
    
    def save_data(self):
        """
        Save the data into .csv files
        """
        self.stns = []
        self.Eframes = []
        self.Bframes = {}
        self.Vframes = []
        for b in self.bins:
            bd = vars(self.bin_details)[str(b)]
            cs = bd.cable_sec
            cs.E.to_csv(self.bdir + "iEt-Bin%02d.csv"%b, float_format="%g")
            cs.B.to_csv(self.bdir + "Bt.csv", float_format="%g")
            cs.V.to_csv(self.bdir + "iVt-Bin%02d.csv"%b, float_format="%g")
            self.Bframes[bd.stn] = cs.B
            self.Eframes.append(cs.E)
            self.stns.append(bd.stn)
            self.Vframes.append(cs.V)
        self.tot_params.to_csv(self.bdir + "sim-params.csv", float_format="%g")
        if self.plot: self.plot_data()
        return
    
    def plot_data(self):
        """
        Plot the corresponding figures.
        """
        base = self.bdir + "/plots/"
        os.makedirs(base, exist_ok=True)
        plotlib.plot_Bxy_stack(list(set(self.stns)), self.Bframes, fbase=base)
        plotlib.plot_Exy_stack(self.stns, self.Eframes, fbase=base)
        plotlib.plot_induced_potential(self.stns, self.Vframes, fbase=base)
        plotlib.plot_total_potential(self.tot_params, fbase=base)
        for b in self.bins:
            fname = self.bdir + "/plots/BExy.Field.B%d.png"%b
            bd = vars(self.bin_details)[str(b)]
            cs = bd.cable_sec
            plotlib.plot_BExy(bd.stn, cs.B, cs.E, fname)
        return