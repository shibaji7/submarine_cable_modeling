"""
    simulate_synB_synT.py: Module is used to implement Synthetic B-Field structures.
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
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from loguru import logger

from gic.model.synthetic import CreateDataSet
from gic.model.oml import OceanModel
from gic.model.cables import CableSection, Cable
from gic.model.conductivity import fetch_static_conductivity_profiles
from gic.model import utils
from gic.model.plotlib import BfieldSummary

class CableAnalysis(object):
    """
    This class is dedicated to B-field simulations and
    storing the results. Class assumes one cable section for a cable
    containing elecetrically long or short cable / cable section.
    This simulation considers
    1. Read magnetic field using
       B-field data_source information.
    2. For each cable segment estimate physical and
       electrical properties, Transfer function,
       induced electric field, and voltage.
    """
    
    def __init__(self, Bfield, cable, out_dir, components=None):
        """
        Parameter:
        ----------
        Bfield: Details of B-field analysis
        p: Tapering value (0-1)
        cable: Cable details holding cable sections
        components: E-field components for nodal analysis and V-calculation
        """
        os.makedirs(out_dir, exist_ok=True)
        self.Bfield = Bfield
        self.cable = cable
        self.p = Bfield.tapering
        self.components = components
        self.out_dir = out_dir
        logger.info(f"B-field run parameters")
        return
    
    def run(self):
        """
        Run all the steps
        """
        self.ocean_model_list =[]
        self.syn = False
        # Genarting / reading all B-field dataset
        self.ds = CreateDataSet(data_sources=self.Bfield.data_sources, dtype="B", p=self.p)
        bdir = self.out_dir
        self.components = self.ds.components if self.components is None else self.components
        self.cbl = Cable(
            cable=self.cable, 
            Efields=None, 
            Bfields=self.ds.field, 
            components=self.components, 
            out_dir=self.out_dir
        )
        self.cbl.setup()
        self.cbl.run_nodal_analysis()
        return
    
    

class SytheticCableAnalysis(object):
    """
    This class is dedicated to synthetic B-field simulations and
    storing the results. Class assumes one cable section for a cable
    containing elecetrically long or short cable / cable section.
    This simulation considers
    1. Create a synthetic magnetic field using
       B-field structure information.
    2. For each cable segment estimate physical and
       electrical properties, Transfer function,
       induced electric field, and voltage.
    """

    def __init__(self, Bfield, cable, out_dir):
        """
        Parameter:
        ----------
        Bfield: Details of B-field analysis
        p: Tapering value (0-1)
        cable: Cable details holding cable sections
        """
        self.Bfield = Bfield
        self.cable = cable
        self.p = Bfield.tapering
        self.out_dir = out_dir
        logger.info(f"Synthetic B-field run parameters")
        return

    def run(self):
        """
        Run all the steps
        """
        self.ocean_model_list =[]
        self.syn = False
        # Genarting / reading all B-field dataset
        if hasattr(self.Bfield, "structure"):
            self.syn = True
            structure = self.Bfield.structure
            self.Am, self.Tm, self.Phim, self.t = (
                np.array(structure.Am),
                np.array(structure.Tm_min) * 60,
                np.array(structure.Phim),
                int(np.array(structure.T_hours) * 60 * 60),
            )
            structure.Am, structure.Tm_min, structure.Phim, structure.t = (
                self.Am,
                self.Tm,
                self.Phim,
                self.t,
            )
            self.ds = CreateDataSet(structure=structure, dtype="B", p=self.p)
        bdir = self.out_dir
        for i, cable_section in enumerate(self.cable.cable_sections):
            cs_physical_prop = CableSection(cable_section.sec_id, cable_section.len_km)
            setattr(self.cable.cable_sections[i], "cs_physical_prop", cs_physical_prop)
            self.Efield = {}
            if hasattr(cable_section.elec_params, "earth_model"):
                earth_model, ocean_depth, ocean_resistivity = (
                    cable_section.elec_params.earth_model,
                    cable_section.elec_params.ocean_depth,
                    cable_section.elec_params.ocean_resistivity,
                )
                earth_profile = fetch_static_conductivity_profiles(earth_model)
                ocean_model = {"depth": ocean_depth, "rho": ocean_resistivity}
                om = OceanModel(
                    thikness=earth_profile["thickness"],
                    resistivity=earth_profile["resistivity"],
                    model_name=earth_model,
                    ocean_model=ocean_model,
                )
                self.ocean_model_list.append(om)
                logger.info(f"Synthetic B {earth_model}->OM({om.model_name})")
                setattr(self.cable.cable_sections[i], "cs_oml", om)
                prep = "h-%d.r-%.2f" % (ocean_depth, ocean_resistivity)
                for stn in self.ds.stations:
                    Enum = self.solve_numerical_Et(om, stn)
                    Km, Tf = self.draw_Km_table(om), self.draw_TF_table(om)
                    Eanl = self.solve_analytical_Et(Tf, stn)
                    setattr(self.cable.cable_sections[i], "cs_Eanl_%s"%stn, Eanl)
                    setattr(self.cable.cable_sections[i], "cs_Enum_%s"%stn, Enum)
                    r = self.check_analytical_numerical(Eanl, Enum)
                    setattr(self.cable.cable_sections[i], "cs_r_%s"%stn, r)
                    pname = bdir + "summary_plot_%s_%s_%s.png" % (
                        earth_model,
                        prep,
                        stn,
                    )
                    self.summary_plot(om, Eanl, Enum, pname, stn)
                    Km.to_csv(
                        bdir + "Kfm_%s_%s_%s.csv" % (earth_model, prep, stn),
                        float_format="%g",
                        header=True,
                        index=False,
                    )
                    Tf.to_csv(
                        bdir + "TF_%s_%s_%s.csv" % (earth_model, prep, stn),
                        float_format="%g",
                        header=True,
                        index=False,
                    )
                    Eanl.to_csv(
                        bdir + "Eanl_%s_%s_%s.csv" % (earth_model, prep, stn),
                        float_format="%g",
                        header=True,
                        index=False,
                    )
                    np.savetxt(
                        bdir + "cor_%s_%s_%s.csv" % (earth_model, prep, stn),
                        np.array([r]).T,
                        delimiter=",",
                        header="rho",
                        comments="",
                    )
                    Enum.to_csv(
                        bdir + "Enum_%s_%s_%s.csv" % (earth_model, prep, stn),
                        float_format="%g",
                        header=True,
                        index=False,
                    )
                    self.ds.field[stn].to_csv(
                        bdir + "Bt_%s.csv" % stn,
                        float_format="%g",
                        header=True,
                        index=False,
                    )
        return

    def draw_Km_table(self, qx):
        """
        Create the Km(f) tables at specific frequencies.
        """
        fm = 1.0 / self.Tm
        Kfm = qx.calcZ(layer="floor", freqs=fm)
        o = pd.DataFrame()
        o["fm (Hz)"], o["Amplitude |Km| (mV/km/nT)"], o["Phase (deg)"] = (
            fm,
            np.absolute(Kfm),
            np.angle(Kfm, deg=True),
        )
        o.index.name = "m"
        return o.copy()

    def draw_TF_table(self, qx):
        """
        Create the TF(f) tables at specific frequencies.
        """
        fm = 1.0 / self.Tm
        o = qx.get_TFs(freqs=fm)
        o["Amplitude (mV/km/nT)"], o["Phase (deg)"] = np.absolute(o.E2B), np.angle(
            o.E2B, deg=True
        )
        o = o.rename(columns={"freq": "fm (Hz)"}).drop(columns=["E2B"])
        o.index.name = "m"
        return o.copy()

    def solve_analytical_Et(self, tf, stn):
        """
        Solve for analytical Et from Bt and TF (tf)
        """
        Bfield = self.ds.field
        Et = np.zeros(len(Bfield[stn]["X"]))
        m = 0
        t = np.array(Bfield[stn].dTime)
        for A, Phi, T in zip(self.Am, self.Phim, self.Tm):
            Et += (
                np.absolute(tf["Amplitude (mV/km/nT)"])[m]
                * A
                * np.sin(2 * np.pi * t / T + np.deg2rad(Phi + tf["Phase (deg)"][m]))
            )
            m += 1
        o = pd.DataFrame()
        o["dTime"], o["X"] = Bfield[stn].dTime, Et
        return o.copy()

    def solve_numerical_Et(self, om, stn):
        """
        Solve for Et from FFT and IFFT operations
        """
        Bfield = self.ds.field
        o = pd.DataFrame()
        o["dTime"] = Bfield[stn].dTime
        for c in self.ds.components:
            Bt, dT = (
                np.array(Bfield[stn][c]),
                Bfield[stn].dTime.iloc[1] - Bfield[stn].dTime.iloc[0],
            )
            Bf, f = utils.fft(Bt, dT)
            E2B = np.array(om.get_TFs(freqs=f).E2B)
            Et = utils.ifft(E2B * Bf)
            o[c] = Et
        return o.copy()

    def check_analytical_numerical(self, Eanl, Enum):
        """
        Regression check for numerical and
        """
        L = int(len(Eanl) / 3)
        for c in self.ds.components:
            r, _ = pearsonr(Eanl[c].tolist()[L:-L], Enum[c].tolist()[L:-L])
            logger.info(f"Corr(Eanl,Enum){': %.10f'%r}")
        return r

    def summary_plot(self, om, Eanl, Enum, pname, stn):
        """
        Create summary plots
        """
        Bfield = self.ds.field
        summ = BfieldSummary()
        ax = summ.add_Bfield_Seq(Bfield[stn], Enum)
        ax.text(
            0.01, 1.05, "Numerical Sol.", ha="left", va="center", transform=ax.transAxes
        )
        ax = summ.add_Bfield_Seq(Bfield[stn], Eanl)
        ax.text(
            0.01,
            1.05,
            "Analytical Sol.",
            ha="left",
            va="center",
            transform=ax.transAxes,
        )
        L = int(len(Eanl) / 3)
        summ.add_Es(Eanl.X.tolist()[L:-L], Enum.X.tolist()[L:-L])
        summ.add_TK_param(om.get_TFs())
        summ.save(pname)
        summ.close()
        return
