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

from OM import OceanModel
import utility

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
        self.eastern_edge = east_edge
        self.western_edge = west_edge
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
        for a in self.components:
            Bt = np.array(self.B[a])
            dT = (self.B.index.tolist()[1]-self.B.index.tolist()[0]).total_seconds()
            Bf, f = utility.fft(Bt, dT)
            E2B = np.array(self.om.get_TFs(freqs=f).E2B)
            Et = utility.ifft(E2B*Bf)
            self.E[a] = Et
        self.E = self.E.set_index("Time")
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
    
    def __init__(self, o, v):
        self.v = v
        self = utility.set_dict(self, o)
        if self.v: 
            logger.info(f"Event analysis and cable parameters")
            utility.print_rec(self)
        ol = pd.read_csv(self.ocean_layers_csv)
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
        return
    
    def calclulate_total_parameters(self):
        """
        Method is responsible for total E-field, V, and 
        V[e-e] caclulations
        """
        return