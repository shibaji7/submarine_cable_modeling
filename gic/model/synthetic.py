"""
    synthetic.py: Module is used to implement Synthetic B-Field structures.
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

# Import required packages
import bezpy
import numpy as np
from scipy import constants as C
import datetime as dt
import pandas as pd
from scipy.stats import pearsonr
from loguru import logger

from .oml import OceanModel
from .plotlib import BfieldSummary


class Preprocess(object):
    """
    This class is used to process the B-field and E-field data.
    """

    def __init__(self, t, field, p=0.1):
        """
        Parameters:
        -----------
        t: Array of time seconds
        field: Array containing the data
        p: Tapering value (0-1)
        """
        self.t = t
        self.field = field
        self.p = p
        return

    def get_tapering_function(self, p=None):
        """
        This method is resposible for generateing
        tapering function based on time sequence t
        and tapering coefficient p
        """
        p = self.p if p is None else p
        T = len(self.t)
        P, P2 = int(T * p), int(T * p / 2)
        w = np.zeros_like(self.t)
        w[:P2] = 0.5 * (1 - np.cos(2 * np.pi * self.t[:P2] / P))
        w[P2 : T - P2] = 1.0
        w[T - P2 :] = 0.5 * (
            1 - np.cos(2 * np.pi * (self.t[-1] - self.t[T - P2 :]) / P)
        )
        return w

    def detrend_magnetic_field(self, p=None):
        """
        This method is resposible for detrend
        magnetic field data and taper it to reduce
        spurious frequency components.
        """
        p = self.p if p is None else p
        fmed = np.median(self.field[:120])
        w = self.get_tapering_function(p)
        f = (self.field - fmed) * w
        return f


class CreateDataSet(object):
    """
    Class is dedicated to generate synthetic dataset for
    magnetic fields. The functionalities includes
    1. Generate synthetic data (B)
    2. Read data
    3. Preprocesss the data
    """

    def __init__(self, data_sources=None, structure=None, dtype="B", p=0.1):
        """
        Parameter:
        ----------
        dataset: Either contains source of dataset or details of generation steps
        """
        self.data_sources = data_sources
        self.structure = structure
        self.p = p
        self.dtype = dtype
        # Compile Dataset
        if data_sources:
            # Read Dataset
            logger.info(f"Reading {self.dtype}-field data")
            if dtype == "B":
                self.read_Bfield_data()
        elif structure:
            # Create Dataset
            logger.info(f"Creating {self.dtype}-field data")
            self.create_synthetic_field(structure)
        # Detranding dataset
        for stn in self.stations:
            for c in self.components:
                proc = Preprocess(
                    np.array(self.field[stn].dTime), np.array(self.field[stn][c])
                )
                self.field[stn][c] = proc.detrend_magnetic_field(p)
        return

    def create_synthetic_field(self, structure):
        """
        This function is responsible for creating a
        synthetic magnetic field with following sets
        of parameters.

        Parameter:
        ----------
        Am (list) - Magntitue at different freuency components (m)
        Phim (list) - Phase at different freuency components (m)
        Tm (list) - Periods of different freuency components (m)
        t (float) - Total time
        """
        self.field = {}
        Am = structure.Am
        Phim = structure.Phim
        Tm = structure.Tm_min
        t = np.linspace(0, structure.t, structure.t, endpoint=False)
        f = np.zeros(len(t))
        for A, Phi, T in zip(Am, Phim, Tm):
            f += A * np.sin(2 * np.pi * t / T + np.deg2rad(Phi))
        self.field["syn"] = pd.DataFrame()
        self.field["syn"]["X"], self.field["syn"]["dTime"] = f, t
        self.components, self.stations = ["X"], ["syn"]
        return

    def read_Bfield_data(self):
        """
        Read B-field dataset
        """
        self.field = {}
        sources = self.data_sources.sources
        self.stations = self.data_sources.stns
        for f in sources:
            stn = f.split("/")[-1].split("_")[0]
            if stn not in self.field.keys():
                self.field[stn] = pd.DataFrame()
            if os.path.exists(f):
                self.field[stn] = pd.concat([self.field[stn], bezpy.mag.read_iaga(f)])
        self.components = ["X", "Y"]
        for stn in self.stations:
            Bf = self.field[stn]
            t0 = Bf.index.tolist()[0].to_pydatetime()
            Bf["dTime"] = (
                np.array(
                    [
                        (x.to_pydatetime() - t0).total_seconds()
                        for x in Bf.reset_index().Time.tolist()
                    ]
                )
            )
            self.field[stn] = Bf
        return
