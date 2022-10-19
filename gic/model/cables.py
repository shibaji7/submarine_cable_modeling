"""cables.py: Module is used to implement cable section analysis and an event study"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import json
import uuid
from types import SimpleNamespace

import bezpy
import numpy as np
import pandas as pd
from loguru import logger

from .conductivity import fetch_static_conductivity_profiles
from .oml import OceanModel
from .plotlib import *
from .utils import *


class TheoriticalSolver(object):
    """
    This class is a theoritcal solution for
    electrically long / short cable-sections.
    """

    def __init__(self, gma, L, E):
        """
        Parameters:
        ------------
        gma: Propagation constant in /km
        L: Cable length in km
        E: Induced electric field in mV/km
        """
        self.gma = gma
        self.L = L
        self.E = E
        return

    def calculate_electrically_short(self, xs=10001):
        """
        Compute voltage along the cable length
        """
        x = np.linspace(0, self.L, xs)  # in km
        U = self.L * self.E / 2  # in Volts
        Ux = U * (2 + self.gma * self.L) * (2 * x - self.L) / (2 * self.L)  # in Volts
        return x, Ux

    def calculate_electrically_long(self, xs=10001):
        """
        Compute voltage along the cable length
        """
        x = np.linspace(0, self.L, xs)  # in km
        U = (1.0 / self.gma) * self.E  # in Volts
        Ux = U * (np.exp(-self.gma * (self.L - x)) - np.exp(-self.gma * x))  # in Volts
        return x, Ux


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

    def __init__(
        self,
        sec_id=None,
        len_km=None,
        directed_length=None,
        elec_params=None,
        loc_i=None,
        loc_f=None,
    ):
        self.sec_id = sec_id if sec_id else uuid.uuid1()
        self.len_km = len_km
        self.directed_length = directed_length
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
        if (
            (self.len_km is None)
            and self.check_location(self.loc_i)
            and self.check_location(self.loc_f)
        ):
            # PARAMETERS OF THE WGS84 EARTH MODEL
            a = 6378.137  # Equatorial radius
            b = 6356.752  # Polar Radius
            e = np.sqrt(0.00669437999014)  # Eccentricity
            lamb = 0.5 * (self.loc_i.lat + self.loc_f.lat)
            self.ln = (111.133 - 0.56 * np.cos(np.deg2rad(2 * lamb))) * np.abs(
                self.loc_f.lat - self.loc_i.lat
            )
            self.le = (
                (111.5065 - 0.1872 * np.cos(np.deg2rad(2 * lamb)))
                * np.cos(np.deg2rad(lamb))
                * np.abs(self.loc_i.lon - self.loc_f.lon)
            )
            self.len_km = np.sqrt(self.ln**2 + self.le**2)
            self.components = ["X", "Y"]
            self.cable_lengths = {"X": self.ln, "Y": self.le}
        elif self.directed_length:
            logger.info("Cable length from directed length")
            self.ln = (
                self.directed_length.ln if hasattr(self.directed_length, "ln") else 0.0
            )
            self.le = (
                self.directed_length.le if hasattr(self.directed_length, "le") else 0.0
            )
            self.components = []
            if hasattr(self.directed_length, "ln"):
                self.components.append("X")
            if hasattr(self.directed_length, "le"):
                self.components.append("Y")
            self.len_km = np.sqrt(self.ln**2 + self.le**2)
            self.cable_lengths = {"X": self.ln, "Y": self.le}
        elif self.len_km:
            self.ln, self.le = self.len_km / np.sqrt(2), self.len_km / np.sqrt(2)
            self.components = ["X", "Y"]
            self.cable_lengths = {"X": self.ln, "Y": self.le}
        else:
            logger.warning("No cable edge information available")
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

    def __init__(
        self,
        sec_id=None,
        len_km=None,
        directed_length=None,
        elec_params=None,
        loc_i=None,
        loc_f=None,
        W=1.0,
        active_termination=False,
    ):
        """
        Properties:
        -----------
        sec_id: Unique cable section ID
        len_km: Length of the cable
        directed_length: Length along a specific direction
        loc_i: Lat, Lon location of the one end of the cable
        loc_f: Lat, Lon location of the another end of the cable
        elec_params: Electrical properties of the cable
        """
        # Compute phyiscal properties of the cable section
        super().__init__(sec_id, len_km, directed_length, elec_params, loc_i, loc_f)
        self.elec_params = elec_params
        self.W = W
        self.active_termination = active_termination
        # Extract electrical properties of the cable
        self.extract_electrical_properties()
        # Extract electrical properties of the cable
        self.compile_oml()
        self.calc_trasmission_line_parameters()
        return

    def compile_oml(self):
        """
        Create ocean model
        """
        if hasattr(self.elec_params, "mtc_model"):
            re = np.array(self.elec_params.mtc_model["resistivity"])
            th = np.array(self.elec_params.mtc_model["thickness"])
            stie = bezpy.mt.Site1d(name="MTC", thicknesses=th[1:], resistivities=re[1:])
            self.om = OceanModel(
                "MTC", None, ocean_model={"depth": th[0], "rho": re[0]}, site=site
            )
        elif (
            hasattr(self.elec_params, "earth_model")
            and hasattr(self.elec_params, "ocean_depth")
            and hasattr(self.elec_params, "ocean_resistivity")
        ):
            d, r, model = (
                self.elec_params.t[0],
                self.elec_params.r[0],
                self.elec_params.earth_model,
            )
            self.om = OceanModel(
                self.elec_params.t[1:],
                self.elec_params.r[1:],
                model_name=model,
                ocean_model={"depth": d, "rho": r},
            )
            logger.info(f"Synthetic {self.sec_id} {model}->OM({self.om.model_name})")
        return

    def extract_electrical_properties(self):
        """
        Extract electrical proerties
        """
        self.electrical_properties = {}
        if hasattr(self.elec_params, "mtc_model"):
            self.electrical_properties["resistivities"] = np.array(
                self.elec_params.mtc_model["resistivity"]
            )
            self.electrical_properties["thicknesses"] = np.array(
                self.elec_params.mtc_model["thickness"]
            )
        else:
            if hasattr(self.elec_params, "r") and hasattr(self.elec_params, "t"):
                self.electrical_properties["resistivities"] = self.elec_params.r
                self.electrical_properties["thicknesses"] = self.elec_params.t
        self.electrical_properties = SimpleNamespace(**self.electrical_properties)
        return

    def calc_trasmission_line_parameters(self):
        """
        Compute the transmission line parameters
        """
        # self.W /= 1e3
        logger.info(f"Cable width: {self.W}")
        if self.elec_params:
            ep = self.electrical_properties
            self.C = self.W * (
                (ep.thicknesses[1] / ep.resistivities[1])
                + (ep.thicknesses[0] / ep.resistivities[0])
            )  # in m/ohm
            self.R = (
                (ep.thicknesses[2] * ep.resistivities[2])
                + (ep.thicknesses[3] * ep.resistivities[3])
            ) / self.W  # in m*ohm
            self.Z, self.Y = 1.0 / self.C, 1.0 / self.R  # in Ohm-m and S/m
            self.gma, self.Z0 = np.sqrt(self.Z * self.Y), np.sqrt(
                self.Z / self.Y
            )  # in /m and Ohm
            if self.active_termination: self.Yn = 1./self.Z0
        else:
            logger.warning("No electrical information available")
        return

    def return_properties(self):
        """
        Create a string of properties for display
        """
        o = "Z: %s (Ohm/m)\n" % (frexp102str(self.Z))
        o += "Y: %s (S/m)\n" % (frexp102str(self.Y))
        o += "Z0: %s (Ohm)\n" % (frexp102str(self.Z0))
        o += "gma: %s (/m)\n" % (frexp102str(self.gma))
        o += "Ad: %s (m)" % (frexp102str(1.0 / self.gma))
        return o

    def compute_eqv_pi_circuit(self, dE, components):
        """
        Calculate equivalent pi circuit model.
        X component is Nort (n), Y component East (e)
        dE: Dataframe containing E-field
        components: [X and Y] for E-fields
        """
        self.Ye, self.Yp2, self.Ie = {}, {}, {}
        if self.active_termination: self.Jn = {}
        for a in components:
            L = self.cable_lengths[a]
            L *= 1000.0  # Convert km to m
            E = np.array(dE[a]) * 1.0e-3 / 1.0e3  # Assuming input mV/km convert to V/m
            self.Ye[a] = 1.0 / (self.Z0 * np.sinh(self.gma * L))
            self.Yp2[a] = (np.cosh(self.gma * L) - 1) * self.Ye[a]
            self.Ie[a] = E / self.Z
            if self.active_termination: self.Jn[a] = E/self.Z
        self.Efield = dE
        self.components = components
        self.compute_Vj(dE.index.tolist())
        return

    def compute_numerical_Et(self, Bfield, components):
        """
        Compute Et using numerical FFT and IFFT block
        Bfield: Dataframe containing B-field
        components: [X and Y] for B-fields
        """
        self.Bfield = Bfield
        self.Efield = pd.DataFrame()
        self.Efield["Time"] = self.Bfield.index.tolist()
        stime = self.Efield.Time.tolist()[0]
        if isinstance(stime, dt.datetime):
            t = np.array(self.Efield.Time.apply(lambda x: (x - stime).total_seconds()))
            self.Efield["dTime"] = t
        else:
            t = np.array(self.Efield.Time)
        for a in components:
            Bt = np.array(self.Bfield[a])
            # Bt = utility.detrend_magnetic_field(np.array(self.Bfield[a]), t)
            dT = t[1] - t[0]
            Bf, f = fft(Bt, dT)
            E2B = np.array(self.om.get_TFs(freqs=f).E2B)
            Et = 2 * ifft(
                component_sign_mappings(
                    "B%sE%s" % (a.lower(), component_mappings("B2E", a).lower())
                )
                * E2B
                * Bf
            )
            self.Efield[component_mappings("B2E", a)] = Et
        self.Efield = self.Efield.set_index("Time")
        # Transforming to E-field components
        self.components = [component_mappings("B2E", a) for a in components]
        self.compute_eqv_pi_circuit(self.Efield, self.components)
        return

    def compute_Vj(self, time):
        """
        Calculate total electric potential induced along the cable segment.
        Vj = Ej_n(t)Lj_n + Ej_e(t)Lj_e
        """
        self.V = pd.DataFrame()
        self.V["Time"] = time
        self.V["Vj"] = 0.0
        for a in self.components:
            lx = self.cable_lengths[a]
            self.V["Vj"] += (
                np.array(self.Efield[a]) * lx
            )  # Potential in mV: E(mV/km) length: km
        self.V = self.V.set_index("Time")
        return

    def calculate_potential_along_cable_section(self, Vi, Vk, ln=1000):
        """
        Caclulate potentials along the cable section
        """
        L = self.len_km
        L *= 1000.0
        x = np.linspace(0, L, ln + 1)
        V = (
            (Vk * np.exp(self.gma * L) - Vi)
            * np.exp(-self.gma * (L - x))
            / (np.exp(self.gma * L) - np.exp(-self.gma * L))
        ) + (
            (Vi * np.exp(self.gma * L) - Vk)
            * np.exp(-self.gma * x)
            / (np.exp(self.gma * L) - np.exp(-self.gma * L))
        )
        return V, x / 1.0e3


class Cable(object):
    """
    This class holds a cable
    Parameters:
    -----------
    cable: Cable parameters
    Efields: Dataframe of E field
    Bfields: Dataframe for B Field
    components: Components for B or E fields
    """

    def __init__(self, cable, Efields, Bfields, components, out_dir):
        self.cable = cable
        self.Efields = Efields
        self.Bfields = Bfields
        self.components = components
        self.out_dir = out_dir
        self.setup()
        return

    def load_electrical_params(self, elec_params):
        """
        Load ocean model electrical parameters
        """
        rf = fetch_static_conductivity_profiles(elec_params.earth_model)
        t, r = np.array(rf["thickness"]), np.array(rf["resistivity"])
        t, r = np.insert(t, 0, elec_params.ocean_depth), np.insert(
            r, 0, elec_params.ocean_resistivity
        )
        setattr(elec_params, "t", t)
        setattr(elec_params, "r", r)
        return elec_params

    def load_LITHO_model(self, loc_i, loc_f, mtc_model):
        cf = ConductivityProfile()
        latlons = np.array([[loc_i.lat, loc_i.lon], [loc_f.lat, loc_f.lon]])
        setattr(elec_params, "mtc_model", cf.compile_bin_profiles(latlons)[0])
        return elec_params

    def setup(self):
        """
        Setup full cable, individual cable segment and
        transmission line
        """
        self.tx_lines = []
        logger.warning("Into cable setup section")
        for j, c in enumerate(self.cable.cable_sections):
            sec_id, len_km, directed_length, elec_params, loc_i, loc_f, stn = (
                None,
                None,
                None,
                None,
                None,
                None,
                "syn",
            )
            active_termination = False
            if hasattr(c, "sec_id"):
                sec_id = c.sec_id
            if hasattr(c, "len_km"):
                len_km = c.len_km
            if hasattr(c, "elec_params"):
                elec_params = c.elec_params
            if hasattr(c, "edge_loc"):
                loc_i, loc_f = c.edge_loc.ini, c.edge_loc.fin
            if hasattr(c, "station"):
                stn = c.station
            if hasattr(c, "directed_length"):
                directed_length = c.directed_length
            if elec_params and hasattr(elec_params, "earth_model"):
                elec_params = self.load_electrical_params(elec_params)
            if elec_params and hasattr(elec_params, "mtc_model"):
                elec_params = self.load_LITHO_model(loc_i, loc_f, elec_params)
            if elec_params and hasattr(elec_params, "active_termination"):
                active_termination = elec_params.active_termination
            tl = TransmissionLine(
                sec_id,
                len_km,
                directed_length,
                elec_params,
                loc_i,
                loc_f,
                1.0,
                active_termination = active_termination
            )
            if self.Efields is not None:
                tl.compute_eqv_pi_circuit(self.Efields[stn], self.components)
                self.eFieldComponents = self.components
            elif self.Bfields is not None:
                tl.compute_numerical_Et(self.Bfields[stn], self.components)
                # Converting to E field components
                self.eFieldComponents = [
                    component_mappings("B2E", a) for a in self.components
                ]
            self.tx_lines.append(tl)
        return

    def run_nodal_analysis(self):
        """
        Run nodal analysis for the cable
        """
        self.nodal_analysis = NodalAnalysis(self.tx_lines, self.eFieldComponents)
        self.nodal_analysis.equivalent_nodel_analysis()
        self.nodal_analysis.solve_admitance_matrix()
        self.noa_result = self.nodal_analysis.consolidate_final_result()
        U0, U1 = None, None
        for a in self.eFieldComponents:
            if U0 is None:
                U0, U1 = self.nodal_analysis.get_voltage_ends_of_cable(a)
            else:
                u0, u1 = self.nodal_analysis.get_voltage_ends_of_cable(a)
                U0 += u0
                U1 += u1
        # Total parameter calculations
        self.tot_params = pd.DataFrame()
        self.tot_params["Time"] = self.tx_lines[0].Efield.index.tolist()
        self.tot_params["V(v)"] = 0.0
        for a in self.eFieldComponents:
            self.tot_params["E." + a] = 0.0
            for i, tl in enumerate(self.tx_lines):
                self.tot_params["E.%s.%02d" % (a, i)] = np.array(tl.Efield[a])
                self.tot_params["E." + a] += np.array(tl.Efield[a])
        for i, tl in enumerate(self.tx_lines):
            self.tot_params["V(v).%02d" % (i)] = np.array(tl.V.Vj) / 1e3
            self.tot_params["V(v)"] += np.array(tl.V.Vj) / 1e3

        self.tot_params["Vt(v)"] = U0 - U1 + np.array(self.tot_params["V(v)"])
        self.tot_params["U0"], self.tot_params["U1"] = U0, U1
        self.tot_params = self.tot_params.set_index("Time")
        self.save_data()
        if ("cable_pot_plot_index" in self.cable.__dict__.keys()) and (
            self.cable.cable_pot_plot_index >= 0
        ):
            self.plots()
        return

    def plots(self):
        bdir = self.out_dir
        Va, La = [], []
        idx = self.cable.cable_pot_plot_index
        for l, tx in enumerate(self.tx_lines):
            pname = bdir + "VCable%02d.png" % l
            U0, U1 = self.nodal_analysis.get_voltage_ends_of_cable_section(
                l, self.eFieldComponents[0]
            )
            if len(self.eFieldComponents) == 2:
                u0, u1 = self.nodal_analysis.get_voltage_ends_of_cable_section(
                    l, self.eFieldComponents[1]
                )
                U0 += u0
                U1 += u1
            U0, U1 = U0[idx], U1[idx]
            V, Lx = tx.calculate_potential_along_cable_section(U0, U1)
            potential_along_section(
                V, Lx, pname, l + 1, U0, U1, tx.Z, tx.Y, tx.gma, tx.Z0
            )
            Va.extend(V.tolist())
            if l == 0:
                La = Lx.tolist()
            else:
                La.extend((Lx + La[-1]).tolist())
        cable_potential(Va, La, bdir + "VCable.png")
        return

    def save_data(self):
        """
        Save all analyzed data including
        Nodal Analysis
        """
        bdir = self.out_dir
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
                    Ji = -1.0 * self.tx_lines[nid].Ie[a]
                    Yii[nid : nid + 2] = np.array(
                        [
                            self.tx_lines[nid].Ye[a] + self.tx_lines[nid].Yp2[a],
                            -self.tx_lines[nid].Ye[a],
                        ]
                    )
                elif nid == self.right_edge:
                    Ji = self.tx_lines[-1].Ie[a]
                    Yii[nid - 1 : nid + 1] = np.array(
                        [
                            -self.tx_lines[-1].Ye[a],
                            self.tx_lines[-1].Yp2[a] + self.tx_lines[-1].Ye[a],
                        ]
                    )
                    if self.tx_lines[-1].active_termination: 
                        Yii[nid] = Yii[nid] + self.tx_lines[-1].Yn
                        Ji = Ji - self.tx_lines[-1].Jn[a]
                else:
                    Ji = self.tx_lines[nid - 1].Ie[a] - self.tx_lines[nid].Ie[a]
                    Yii[nid - 1 : nid + 2] = np.array(
                        [
                            -self.tx_lines[nid - 1].Ye[a],
                            self.tx_lines[nid - 1].Ye[a]
                            + self.tx_lines[nid].Ye[a]
                            + self.tx_lines[nid - 1].Yp2[a]
                            + self.tx_lines[nid].Yp2[a],
                            -self.tx_lines[nid].Ye[a],
                        ]
                    )
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
            self.V[a] = np.matmul(iY, J)
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
            o["cables"][bid] = {
                "R": tx.R,
                "C": tx.C,
                "gma": tx.gma,
                "Z0": tx.Z0,
                "ln": tx.ln,
                "le": tx.le,
                "len_km": tx.len_km,
                "Ye": {},
                "Yp2": {},
                "Ie": {},
            }
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
        u = 1.0 if unit == "V" else 1000.0
        U0, U1 = np.round(self.V[comp][0, :] * u, 2), np.round(
            self.V[comp][-1, :] * u, 2
        )
        logger.info(f"Max(V) at the end (Component-{comp}), {np.max(U0)} {np.max(U1)}")
        return U0, U1

    def get_voltage_ends_of_cable_section(self, b=0, comp="X", unit="V"):
        """
        Provide the voltage at the ends of the
        cable to calculate total voltage
        """
        u = 1.0 if unit == "V" else 1000.0
        U0, U1 = np.round(self.V[comp][b, :] * u, 2), np.round(
            self.V[comp][b + 1, :] * u, 2
        )
        logger.info(
            f"Max(V) at the end of Section-{b}(Component-{comp}), {np.max(U0)} {np.max(U1)}"
        )
        return U0, U1


class Node(object):
    """
    Blank node for nodal analysis
    """

    def __init__(self):
        return
