"""DSTL.py: Module is used to implement Distributed Source Transmission Line theory"""

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
import numpy as np
import pandas as pd
import json

from loguru import logger

class TransmissionLine(object):
    """
    This class is dedicated to calculate the Earth's 
    potential across the cable segment.
    """
    
    def __init__(self, model_num, site, bin_details, ocean_model={"depth":5e3, "rho":0.25}, w=1.):
        logger.info(f"Compile TLM[{model_num}] to calc Earth's pot.")
        self.model_num = model_num
        self.ocean_model = ocean_model
        self.bd = bin_details
        self.site = site
        self.w = w
        self.components = ["X", "Y"]
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
        phi = 0.5*(self.bd.western_edge[1] + self.bd.eastern_edge[1])
        Ln = (111.133-0.56*np.cos(np.deg2rad(2*phi)))*np.abs(self.bd.eastern_edge[0]-self.bd.western_edge[0])
        Le = (111.5065-0.1872*np.cos(np.deg2rad(2*phi)))*np.cos(np.deg2rad(phi))*\
                np.abs(self.bd.eastern_edge[1]-self.bd.western_edge[1])
        return Ln, Le
    
    def calc_tx_parameters(self):
        """
        Calculate transmission line parameters:
        Propagation constant (gma), characteristics
        impedance (Z0).
        """
        self.C = self.w*( (self.site.thicknesses[0]/self.site.resistivities[0]) + 
                         (self.ocean_model["depth"]/self.ocean_model["rho"]) )
        self.R = ( (self.site.thicknesses[1]*self.site.resistivities[1]) + 
                  (self.site.thicknesses[2]*self.site.resistivities[2]) )/self.w
        self.Z, self.Y = 1./self.C, 1./self.R
        self.gma, self.Z0 = np.sqrt(self.Z*self.Y), np.sqrt(self.Z/self.Y)
        return
    
    def compute_eqv_pi_circ_model(self, dE):
        """
        Calculate equivalent pi circuit model.
        X component is Nort (n), Y component East (e)
        """
        self.Ln, self.Le = self.segmentL()
        self.Ye, self.Yp2, self.Ie = {}, {}, {}
        for a, L in zip(self.components, [self.Ln, self.Le]):
            E = np.array(dE[a])
            self.Ye[a] = 1./( self.Z0*np.sinh(self.gma*L) )
            self.Yp2[a] = ( np.cosh(self.gma*L)-1 ) * self.Ye[a]
            self.Ie[a] = E/self.Z
        return
    
    def compiles(self, dE):
        """
        Compiles all the small section of computation
        cronologicaly.
        """
        self.calc_tx_parameters()
        self.compute_eqv_pi_circ_model(dE)
        return

class NodalAnalysis(object):
    """
    Nodal analysis of transmission line model.
    Assumption: Pi-models are connecetd sequentially and
    linearly, thus with N cable sections there are 
    N+1 nodes to be analyzed.
    """
    
    def __init__(self, tx_lines, bdir):
        self.bdir = bdir
        logger.info(f"In nodal analysis for {len(tx_lines)+1} nodes")
        self.tx_lines = tx_lines
        self.nnodes = len(tx_lines) + 1
        self.node_ids = np.arange(self.nnodes)
        self.left_edge, self.right_edge = 0, self.node_ids[-1]
        self.nodes = {}
        self.components = ["X", "Y"]
        return
    
    def equivalent_nodel_analysis(self):
        """
        Nodal analysis of the network
        """
        logger.info(f"Eq. nodal analysis.")
        for nid in self.node_ids:
            self.nodes[nid] = {}
            logger.info(f"Node:{nid}")
            for a in self.components:
                node = Node()
                Yii = np.zeros_like(self.node_ids)
                if nid == self.left_edge: 
                    Ji = -1.*self.tx_lines[nid].Ie[a]
                    Yii[nid:nid+2] = np.array([self.tx_lines[nid].Ye[a]+\
                                               self.tx_lines[nid].Yp2[a],
                                               -self.tx_lines[nid].Ye[a]])
                elif nid == self.right_edge: 
                    Ji = self.tx_lines[-1].Ie[a]
                    Yii[nid-1:nid+1] = np.array([-self.tx_lines[-1].Ye[a],
                                                 self.tx_lines[-1].Yp2[a]+\
                                                 self.tx_lines[-1].Ye[a]])
                else: 
                    Ji = self.tx_lines[nid-1].Ie[a] - self.tx_lines[nid].Ie[a]
                    Yii[nid-1:nid+2] = np.array([-self.tx_lines[nid-1].Ye[a], 
                                                 self.tx_lines[nid-1].Ye[a]+\
                                                 self.tx_lines[nid].Ye[a]+\
                                                 self.tx_lines[nid-1].Yp2[a]+\
                                                 self.tx_lines[nid].Yp2[a],
                                                 -self.tx_lines[nid].Ye[a]])
                setattr(node, "Ji", Ji)
                setattr(node, "Yii", Yii)
                self.nodes[nid][a] = node
        return
    
    def solve_admitance_matrix(self):
        """
        Solve: [V] = inv([Y]).[J]
        """
        self.V = {}
        logger.info(f"Solving admitance matrix.")
        for a in self.components:
            logger.info(f"Solving for component {a}.")
            J, Y = [], []
            for nid in self.node_ids:
                n = self.nodes[nid][a]
                J.append(n.Ji)
                Y.append(n.Yii)
            J, Y = np.array(J), np.array(Y)
            logger.info(f"Sh(J):{J.shape}, Sh(Y):{Y.shape}")
            iY = np.linalg.inv(Y)
            self.V[a] = np.matmul(iY,J)
            logger.info(f"Sh(V):{self.V[a].shape}")
        return
    
    def consolidate_final_result(self):
        """
        Estimated Vi,k are the voltages at the end 
        of the cable sections. Here we store estimated 
        voltages in csv format, and line parameters 
        (R, C, gma, L, Z0, Ln, Le, Ye[n,e], Yp2[n,e], Ie[n,e])
        in json format.
        """
        o = {"nodes": {}, "cables": {}}
        logger.info(f"Consolidate all results.")
        for bid, tx in enumerate(self.tx_lines):
            bid += 1
            o["cables"][bid] = {"R": tx.R, "C": tx.C, "gma": tx.gma, "Z0": tx.Z0, "Ln": tx.Ln, 
                                "Le": tx.Le, "Ye": {}, "Yp2": {}, "Ie": {}}
            for a in self.components:
                o["cables"][bid]["Ye"][a] = tx.Ye[a]
                o["cables"][bid]["Yp2"][a] = tx.Yp2[a]
                o["cables"][bid]["Ie"][a] = tx.Ie[a].tolist()
        for nid in self.node_ids:
            nid = str(nid)
            for a in self.components:
                n = self.nodes[int(nid)][a]
                o["nodes"][nid] = {a: {}}
                o["nodes"][nid][a]["Ji"] = n.Ji.tolist()
                o["nodes"][nid][a]["Yii"] = n.Yii.tolist()
        with open(self.bdir + "est_cable_props.json", "w") as f: 
            f.write(json.dumps(o, sort_keys=True, indent=4))
        return
    
    def get_voltage_ends_of_cable(self, comp="X", unit="V"):
        """
        Provide the voltage at the ends of the 
        cable to calculate total voltage
        """
        u = 1.e-3 if unit == "V" else 1.
        U0, U1 = np.round(self.V[comp][0,:]*u,2), np.round(self.V[comp][-1,:]*u,2)
        logger.info(f"Max(V) at the end (Component-{comp}), {np.max(U0)} {np.max(U1)}")
        return U0, U1
    
    
class Node(object):
    """
    Blank node for nodal analysis
    """
    def __init__(self):
        return