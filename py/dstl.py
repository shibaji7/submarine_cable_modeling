"""DSTL.py: Module is used to implement Distributed Source Transmission Line theory"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
os.path.extend(["py/", "py/config/"])
import numpy as np
import pandas as pd

from loguru import logger

class TransmissionLine(object):
    """
    This class is dedicated to calculate the Earth's 
    potential across the cable segment.
    """
    
    def __init__(self, site, cd, ocean_model={"depth":5e3, "rho":0.25}, w=1.):
        logger.info(f"Compile TLM[{model_num}] to calc Earth's pot.")
        self.model_num = model_num
        self.ocean_model = ocean_model
        self.cable_details = cd
        self.site = site
        self.w = w
        return
    
    def calc_tx_parameters(self):
        """
        Calculate transmission line parameters:
        Propagation constant (gma), characteristics
        impedance (Z0).
        """
        C = self.w*( (self.site.thicknesses[0]/self.site.resistivities[0]) + 
                    (self.ocean_model["depth"]/self.ocean_model["rho"]) )
        R = ( (self.site.thicknesses[1]*self.site.resistivities[1]) + 
             (self.site.thicknesses[2]*self.site.resistivities[2]) )/self.w
        self.Z, self.Y = 1./C, 1./R
        self.gma, self.Z0 = np.sqrt(Z*Y), np.sqrt(Z/Y)
        return
    
    def compute_eqv_pi_circ_model(self, E):
        """
        Calculate equivalent pi circuit model
        """
        L = 0.
        self.Ye = 1./( self.Z0*np.sinh(self.gma*L) )
        self.Yp2 = ( np.cosh(self.gma*L)-1 ) * self.Ye
        self.Ie = E/self.Z
        return
    
    def compiles(self):
        """
        Compiles all the small section of computation
        cronologicaly.
        """
        self.calc_tx_parameters()
        self.compute_eqv_pi_circ_model()
        return