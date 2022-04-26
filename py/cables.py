"""cables.py: Module is used to implement cable section analysis and an event study"""

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
import uuid
import numpy as np
from loguru import logger
from types import SimpleNamespace
import bezpy
import pandas as pd
import datetime as dt
import json

import utility
from oml import OceanModel
import plotlib

class ElectricalProps(object):
    """
    Store electrical properties
    """
    
    def __init__(self, path=None, model=None, ):
        return

class CableSection(object):
    """
    This class holds a cable section of a big cable
    Parameters:
    -----------
    sec_id: Cable Section ID
    edge_locations: Edge informations (lat, lon); list of tuple
    electrical_properties: Electrical properties of the cable
    c_len: Length of the cable
    le: Eastward length of the cable
    ln: Northward length of the cable
    """
    
    def __init__(self, sec_id=None, len_km=None, elec_params=None, loc_i=None, loc_f=None):
        self.sec_id = sec_id if sec_id else uuid.uuid1()
        self.len_km = len_km
        self.elec_params = elec_params
        self.loc_i = loc_i
        self.loc_f = loc_f
        self.compute_lengths()
        return
    
    def check_location(self, loc):
        """
        Check lat/lon exists in file or not
        """
        tag = True if (hasattr(loc, "lat") and hasattr(loc, "lon")) else False
        return tag
    
    def compute_lengths(self):
        """
        Compute all the length of the Cable Section
        """
        if (self.len_km is None) and self.check_location(self.loc_i) and self.check_location(self.loc_f):
            # PARAMETERS OF THE WGS84 EARTH MODEL
            a = 6378.137 # Equatorial radius
            b = 6356.752 # Polar Radius
            e = np.sqrt(0.00669437999014) # Eccentricity
            phi = 0.5*(self.loc_i.lon + self.loc_f.lon)
            self.ln = (111.133-0.56*np.cos(np.deg2rad(2*phi)))*\
                np.abs(self.loc_f.lat-self.loc_i.lat)
            self.le = (111.5065-0.1872*np.cos(np.deg2rad(2*phi)))*np.cos(np.deg2rad(phi))*\
                np.abs(self.loc_i.lon-self.loc_f.lon)
            self.len_km = np.sqrt(self.ln**2+self.le**2)
        elif self.len_km: self.ln, self.le = self.len_km/np.sqrt(2), self.len_km/np.sqrt(2)
        else: logger.warning("No cable edge information available")
        return
    
class TransmissionLine(CableSection):
    """
    This class is dedicated for DSTL.
    Parameters:
    -----------
    electrical_properties: Electrical properties of the ground/ocean layers
    W: Reset parameter (1)
    ** All properties/methods of the Cable Sections
    """
    
    def __init__(self, sec_id=None, len_km=None, elec_params=None, 
                 loc_i=None, loc_f=None, W=1., model_path=None):
        super().__init__(sec_id, len_km, elec_params, loc_i, loc_f)
        self.elec_params = elec_params
        self.model_path = model_path
        self.W = W
        self.compile_oml()
        self.extract_electrical_properties()
        self.calc_trasmission_line_parameters()
        return
    
    def compile_oml(self):
        """
        Create ocean model
        """
        if hasattr(self.elec_params, "earth_model")\
            and hasattr(self.elec_params, "ocean_depth")\
            and hasattr(self.elec_params, "ocean_resistivity"):
            d, r, m = self.elec_params.ocean_depth, self.elec_params.ocean_resistivity,\
                    self.elec_params.earth_model
            mpath = self.model_path%"BME" if "uniform" in self.model_path else\
                self.model_path%self.elec_params.earth_model
            self.om = OceanModel(m, mpath, ocean_model={"depth":d, "rho":r})
            if "uniform" in m:
                rho = float(m.split(".")[-1])
                self.om.site.resistivities[:] = rho
            logger.info(f"Synthetic {self.sec_id} {m}->OM({self.om.model_name})")
        return
    
    def extract_electrical_properties(self):
        """
        Extract electrical proerties
        """
        self.electrical_properties = {}
        if hasattr(self.elec_params, "r"): self.electrical_properties["resistivities"] = self.elec_params.r
        if hasattr(self.elec_params, "t"): self.electrical_properties["thicknesses"] = self.elec_params.t
        if hasattr(self.elec_params, "earth_model"):
            mpath = self.model_path%"BME" if "uniform" in self.model_path else\
                self.model_path%self.elec_params.earth_model
            r, t = [], []
            if hasattr(self.elec_params, "ocean_depth"): t.append(self.elec_params.ocean_depth)
            if hasattr(self.elec_params, "ocean_resistivity"): r.append(self.elec_params.ocean_resistivity)
            site = bezpy.mt.read_1d_usgs_profile(mpath)
            if "uniform" in self.elec_params.earth_model:
                rho = float(self.elec_params.earth_model.split(".")[-1])
                site.resistivities[:] = rho
            r.extend(site.resistivities.tolist())
            t.extend(site.thicknesses.tolist())            
            self.electrical_properties["resistivities"] = r
            self.electrical_properties["thicknesses"] = t
        self.electrical_properties = SimpleNamespace(**self.electrical_properties)
        return
    
    def calc_trasmission_line_parameters(self):
        """
        Compute the transmission line parameters
        """
        if self.elec_params:
            ep = self.electrical_properties
            self.C = self.W*( (ep.thicknesses[1]/ep.resistivities[1]) + 
                         (ep.thicknesses[0]/ep.resistivities[0]) )
            self.R = ( (ep.thicknesses[2]*ep.resistivities[2]) + 
                  (ep.thicknesses[3]*ep.resistivities[3]) )/self.W
            self.Z, self.Y = 1./self.C, 1./self.R
            self.gma, self.Z0 = np.sqrt(self.Z*self.Y), np.sqrt(self.Z/self.Y)
        else: logger.warning("No electrical information available")
        return
    
    def compute_eqv_pi_circuit(self, dE, components):
        """
        Calculate equivalent pi circuit model.
        X component is Nort (n), Y component East (e)
        dE: Dataframe containing E-field
        components: [X and Y]
        """
        self.Ye, self.Yp2, self.Ie = {}, {}, {}
        for a, L in zip(components, [self.ln, self.le]):
            L *= 1000. # Convert km to m
            E = np.array(dE[a])
            self.Ye[a] = 1./( self.Z0*np.sinh(self.gma*L) )
            self.Yp2[a] = ( np.cosh(self.gma*L)-1 ) * self.Ye[a]
            self.Ie[a] = E/self.Z
        self.Efield = dE
        self.components = components
        self.compute_Vj(dE.index.tolist())
        return
    
    def compute_numerical_Et(self, Bfield, components):
        """
        Compute Et using numerical FFT and IFFT block 
        """
        self.components = components
        self.Bfield = Bfield
        self.Efield = pd.DataFrame()
        self.Efield["Time"] = self.Bfield.index.tolist()
        stime = self.Efield.Time.tolist()[0]
        if isinstance(stime, dt.datetime):
            t = np.array(self.Efield.Time.apply(lambda x: (x-stime).total_seconds()))
            self.Efield["dTime"] = t
        else: t = np.array(self.Efield.Time)
        for a in self.components:
            Bt = utility.detrend_magnetic_field(np.array(self.Bfield[a]), t)
            dT = (t[1]-t[0])
            Bf, f = utility.fft(Bt, dT)
            E2B = np.array(self.om.get_TFs(freqs=f).E2B)
            Et = 2*utility.ifft(E2B*Bf)
            self.Efield[a] = Et
        self.Efield = self.Efield.set_index("Time")
        self.compute_eqv_pi_circuit(self.Efield, components)
        self.compute_Vj(self.Bfield.index.tolist())
        return
    
    def compute_Vj(self, time):
        """
        Calculate total electric potential induced along the cable segment.
        Vj = Ej_n(t)Lj_n + Ej_e(t)Lj_e
        """
        Ln, Le = self.ln, self.le
        self.V = pd.DataFrame()
        self.V["Time"] = time
        self.V["Vj"] = 0.
        for a in self.components:
            lx = Ln if a == "X" else Le
            self.V["Vj"] = np.array(self.Efield[a])*lx# + np.array(self.Efield.Y)*Le)
        self.V = self.V.set_index("Time")
        return
    
    def calculate_potential_along_cable_section(self, Vi, Vk, ln=1000):
        """
        Caclulate potentials along the cable section
        """
        L = self.len_km
        L *= 1000.
        x = np.linspace(0, L, ln+1)
        V = ( (Vk*np.exp(self.gma*L)-Vi)*np.exp(-self.gma*(L-x))/(np.exp(self.gma*L)-np.exp(-self.gma*L)) ) +\
            ( (Vi*np.exp(self.gma*L)-Vk)*np.exp(-self.gma*x)/(np.exp(self.gma*L)-np.exp(-self.gma*L)) )
        return V, x/1.e3

class Cable(object):
    """
    This class holds a cable
    Parameters:
    -----------
    """
    
    def __init__(self, args, Efields, Bfields, components):
        self.args = args
        self.Efields = Efields
        self.Bfields = Bfields
        self.components = components
        self.setup()
        return
    
    def load_electrical_params(self, elec_params, j):
        """
        Load ocean model electrical parameters
        """
        o = pd.read_csv(self.args.cable.ocean_layers_csv)
        o = o[o["bin"]=="Bin%d"%(j+1)]
        setattr(elec_params, "ocean_depth", o["depth"].tolist()[0])
        setattr(elec_params, "ocean_resistivity", o["rho"].tolist()[0])
        return elec_params
    
    def setup(self):
        """
        Setup full cable, individual cable segment and
        transmission line
        """
        self.tx_lines = []
        logger.warning("Into cable setup section")
        for j, c in enumerate(self.args.cable.cable_sections):
            sec_id, len_km, elec_params, loc_i, loc_f, stn = None, None, None, None, None, "syn"
            if hasattr(c, "sec_id"): sec_id = c.sec_id
            if hasattr(c, "len_km"): len_km = c.len_km
            if hasattr(c, "elec_params"): elec_params = c.elec_params
            if hasattr(c, "edge_loc"): loc_i, loc_f = c.edge_loc.ini, c.edge_loc.fin
            if hasattr(c, "station"): stn = c.station
            if elec_params and hasattr(elec_params, "load_param"): 
                elec_params = self.load_electrical_params(elec_params, j)
            tl = TransmissionLine(sec_id, len_km, elec_params, loc_i, loc_f, 1., self.args.cable.model_dir_path)
            if self.Efields is not None: tl.compute_eqv_pi_circuit(self.Efields[stn], self.components)
            elif self.Bfields is not None: tl.compute_numerical_Et(self.Bfields[stn], self.components)
            self.tx_lines.append(tl)
        return
    
    def run_nodal_analysis(self):
        """
        Run nodal analysis for the cable
        """
        self.nodal_analysis = NodalAnalysis(self.tx_lines, self.components)
        self.nodal_analysis.equivalent_nodel_analysis()
        self.nodal_analysis.solve_admitance_matrix()
        self.noa_result = self.nodal_analysis.consolidate_final_result()
        U0, U1 = None, None
        for a in self.components:
            if U0 is None: U0, U1 = self.nodal_analysis.get_voltage_ends_of_cable(a)
            else: 
                u0, u1 = self.nodal_analysis.get_voltage_ends_of_cable(a)
                U0 += u0; U1 += u1
        # Total parameter calculations
        self.tot_params = pd.DataFrame()
        self.tot_params["Time"] = self.tx_lines[0].Efield.index.tolist()
        self.tot_params["V(v)"] = 0.
        for a in self.components:
            self.tot_params["E."+a] = 0.
            for i, tl in enumerate(self.tx_lines):
                self.tot_params["E."+a] += tl.Efield[a]
        for i, tl in enumerate(self.tx_lines):
            self.tot_params["V(v)"] += tl.V.Vj

        self.tot_params["V(v)"] /= 1000.
        self.tot_params["Vt(v)"] = (U0-U1+np.array(self.tot_params["V(v)"]))/1000.
        self.tot_params = self.tot_params.set_index("Time")
        self.save_data()
        if ("cable_pot_plot_index" in self.args.cable.__dict__.keys()) and\
            (self.args.cable.cable_pot_plot_index >= 0): self.plots()
        return
    
    def plots(self):
        bdir = self.args.out_dir
        Va, La = [], []
        idx = self.args.cable.cable_pot_plot_index
        for l, tx in enumerate(self.tx_lines):
            pname = bdir + "VCable%02d.png"%l
            U0, U1 = self.nodal_analysis.get_voltage_ends_of_cable_section(l, self.components[0])
            if len(self.components) == 2:
                u0, u1 = self.nodal_analysis.get_voltage_ends_of_cable_section(l, self.components[1])
                U0 += u0; U1 += u1
            U0, U1 = U0[idx], U1[idx]
            V, Lx = tx.calculate_potential_along_cable_section(U0, U1)
            plotlib.potential_along_section(V, Lx/1.e3, pname, l+1,
                                            U0, U1, tx.Z, tx.Y, tx.gma, tx.Z0)
            Va.extend(V.tolist())
            if l==0: La = Lx.tolist()
            else: La.extend((Lx+La[-1]).tolist())
        plotlib.cable_potential(Va, La, bdir + "VCable.png")
        return
    
    def save_data(self):
        """
        Save all analyzed data including 
        Nodal Analysis
        """
        bdir = self.args.out_dir
        with open(bdir + "est_cable_props.json", "w") as f: 
            f.write(json.dumps(self.noa_result, sort_keys=True, indent=4))
        self.tot_params.to_csv(bdir + "sim-params.csv", float_format="%g")
        return
    
class NodalAnalysis(object):
    """
    Nodal analysis of transmission line model.
    Assumption: Pi-models are connecetd sequentially and
    linearly, thus with N cable sections there are 
    N+1 nodes to be analyzed.
    """
    
    def __init__(self, tx_lines, components):
        logger.info(f"In nodal analysis for {len(tx_lines)+1} nodes")
        self.tx_lines = tx_lines
        self.nnodes = len(tx_lines) + 1
        self.node_ids = np.arange(self.nnodes)
        self.left_edge, self.right_edge = 0, self.node_ids[-1]
        self.nodes = {}
        self.components = components
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
                Yii = np.zeros_like(self.node_ids, dtype=float)
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
            o["cables"][bid] = {"R": tx.R, "C": tx.C, "gma": tx.gma, "Z0": tx.Z0, "ln": tx.ln, 
                                "le": tx.le, "len_km": tx.len_km, "Ye": {}, "Yp2": {}, "Ie": {}}
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
        return o
    
    def get_voltage_ends_of_cable(self, comp="X", unit="V"):
        """
        Provide the voltage at the ends of the 
        cable to calculate total voltage
        """
        u = 1.e-3 if unit == "V" else 1.
        U0, U1 = np.round(self.V[comp][0,:]*u,2), np.round(self.V[comp][-1,:]*u,2)
        logger.info(f"Max(V) at the end (Component-{comp}), {np.max(U0)} {np.max(U1)}")
        return U0, U1
    
    def get_voltage_ends_of_cable_section(self, b=0, comp="X", unit="V"):
        """
        Provide the voltage at the ends of the 
        cable to calculate total voltage
        """
        u = 1.e-3 if unit == "V" else 1.
        U0, U1 = np.round(self.V[comp][b,:]*u,2), np.round(self.V[comp][b+1,:]*u,2)
        logger.info(f"Max(V) at the end of Section-{b}(Component-{comp}), {np.max(U0)} {np.max(U1)}")
        return U0, U1
    
    
class Node(object):
    """
    Blank node for nodal analysis
    """
    def __init__(self):
        return