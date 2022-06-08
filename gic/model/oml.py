"""
    oml.py: Module is used to implement Ocean Model with Layered-Earth structures.
    Class:
    ------
    OcenModel: Dedicated for ocean - Earth B and E field transfer functions
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import bezpy
import numpy as np
from scipy import constants as C
import pandas as pd
from loguru import logger


class OceanModel(object):
    """
    This class is to emulate electric and magnetic fields
    in a 1D ocean. We provide the 1D Earth conductivity model and
    one layer of ocean parameters (depth and resistivity) to emulate
    1D ocean.
    """

    def __init__(
        self,
        thikness,
        resistivity,
        model_name="",
        ocean_model={"depth": 5e3, "rho": 0.25},
        site=None,
        flim=[1e-6, 1e0],
        kind="floor",
    ):
        """
        Initialize the ocean model
        """
        logger.info(f"Compile OM[{model_name}] to calc O({kind}) E- and B-Fields")
        self.model_name = model_name
        self.ocean_model = ocean_model
        self.kind = kind
        self.freqs = np.linspace(flim[0], flim[1], int(flim[1] / flim[0]) + 1)
        self.functions = {
            "E2B": lambda Z, Zd, kd: Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Z)),
            "B2B": lambda Z, Zd, kd: 1.0 / (np.cosh(kd) + (Zd * np.sinh(kd) / Z)),
        }

        if site:
            self.site = site
        else:
            self.site = bezpy.mt.Site1d(
                name=model_name, thicknesses=thikness, resistivities=resistivity
            )
            if kind != "floor":
                self.site.resistivities = np.insert(
                    self.site.resistivities, 0, ocean_model["rho"]
                )
                self.site.thicknesses = np.insert(
                    self.site.thicknesses, 0, ocean_model["depth"]
                )
        return

    def calcZ(self, layer="ocean", freqs=None):
        """
        Compute the Km(f) for different layers
        """
        freqs = np.copy(self.freqs) if freqs is None else freqs
        omega = 2 * C.pi * freqs
        if self.kind == "floor":
            sigma_s = 1 / self.ocean_model["rho"]
            k2 = 1.0j * omega * C.mu_0 * sigma_s
            k = np.sqrt(k2)
            Zo = (1.0j * omega * C.mu_0 / k) / (C.mu_0 / 1.0e-3)
            Kf = self.site.calcZ(freqs)[1]
            self.Z = {"ocean": Zo, "floor": Kf}
        else:
            self.Z = {
                "ocean": self.site.calcZ(freqs, lev=0)[1],
                "floor": self.site.calcZ(freqs, lev=1)[1],
            }
        return self.Z[layer]

    def calcTF(self, kinds=["E2B", "B2B"], freqs=None):
        """
        Calculate the transfer functions.
        """
        freqs = np.copy(self.freqs) if freqs is None else freqs
        omega = 2 * C.pi * freqs
        Zd = self.calcZ("floor", freqs=freqs)
        Z = self.calcZ("ocean", freqs=freqs)
        TFs = {}

        omega = 2 * C.pi * freqs
        sigma_s = 1 / self.ocean_model["rho"]
        k2 = 1.0j * omega * C.mu_0 * sigma_s
        k = np.sqrt(k2)
        kd = k * self.ocean_model["depth"]

        for kind in kinds:
            TFs[kind] = self.functions[kind](Z, Zd, kd)
        return TFs

    def get_TFs(self, key="E2B", freqs=None):
        """
        Fetch the transfer function based on key
        and frequency.
        """
        TFs = self.calcTF(freqs=freqs)
        tf = pd.DataFrame()
        tf["freq"], tf[key] = np.copy(self.freqs) if freqs is None else freqs, TFs[key]
        return tf
