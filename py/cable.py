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
    
    def __init__(self, model_name, model_path, stn_name, 
                 east_edge, west_edge, data_files):
        self.model_name = model_name
        self.model_path = model_path
        self.stn_name = stn_name
        self.eastern_edge = eastern_edge
        self.western_edge = western_edge
        self.data_files = data_files
        self.load_dataset()
        return
    
    def load_dataset(self):
        """
        Load dataset from given data files
        """
        o = pd.DataFrame()
        for f in self.data_files:
            if os.path.exists(f): o = pd.concat([o, bezpy.mag.read_iaga(f)])
        self.B = o.copy()
        return
    
    def compute_Vj(self):
        """
        Calculate total electric potential induced along the cable segment.
        Vj = Ej_n(t)Lj_n + Ej_e(t)Lj_e
        """
        Ln, Le = self.segmentL()
        self.V = pd.DataFrame()
        self.V["Vj"] = np.array(self.E.X)*Ln + np.array(self.E.Y)*Le
        self.V["Date"] = self.E.index.tolist()
        self.V = self.V.set_index("Date")
        return
    
    def compute_numerical_Et(self):
        """
        Compute Et using numerical FFT and IFFT block 
        """
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
        return