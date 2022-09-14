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

def efmte(x, n=3, m=1, e="e", nmax=16):
    n = 3 if x < 1. else 0
    def _expc(s): # return exponent character of e-formatted float
        return next(filter(lambda character: character in {"E", "e"}, s))
    def _pad0(s, n): # return string padded to length n
        return ("{0:0>" + str(n) + "}").format(s)
    def _efmtes(s, n): # reformat e-formatted float: n e-digits
        m, e, p = s.partition(_expc(s)) # mantissa, exponent, +/-power
        return m + e + p[0] + _pad0(p[1:], n)
    def _efmt(x, n, e): # returns formatted float x: n decimals, "e"/"E"
        return ("{0:." + str(n) + e + "}").format(x)
    x = x if isinstance(x, float) else float("nan")
    nmax = 16 if not isinstance(nmax, int) else max(0, nmax)
    n = 6 if not isinstance(n, int) else min(max(0, n), nmax)
    m = 2 if not isinstance(m, int) else max(0, m)
    e = "e" if e not in {"E", "e"} else e
    return _efmtes(_efmt(x, n, e), m)

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
    
    def plot_earth_interior(self, ax, params={"N": 1001, "t_min":60., 
                                              "t_max":120., "r_min":0., 
                                              "tag_offset":0.4, "thick":2.},
                           fontdict={"size":7}):
        ###########################################
        # Earth interior with Ocean schematics
        ###########################################
        """
        Method is used to plot a schematic of Earth"s schematics with conductivity

        Parameters:
        -----------
        ax: Matplotlib axis object
        """
        import matplotlib.pyplot as plt
        
        N = params["N"] if "N" in params.keys() else 1001
        t_min = params["t_min"] if "t_min" in params.keys() else 60.
        t_max = params["t_max"] if "t_max" in params.keys() else 120.
        r_min =  params["r_min"] if "r_min" in params.keys() else 0.
        tag_offset = params["tag_offset"] if "tag_offset" in params.keys() else 0.4
        thick = params["thick"] if "thick" in params.keys() else 2
        layers = np.arange(1, len(self.site.thicknesses)+2)
        d = thick/len(layers)
        r_max = np.sum(layers)*d/thick
        
        
        ax.plot([0], [0], "w", lw=0)
        ax.grid(False)
        plt.xticks([0], " ")
        plt.yticks([0], " ")
        ax.set_thetamin(t_min)
        ax.set_thetamax(t_max)
        ax.set_rmin(r_min)
        ax.set_rmax(r_max)
        
        resistivities = self.site.resistivities.tolist()
        thicknesses = self.site.thicknesses.tolist()
        resistivities.insert(0, self.ocean_model["rho"])
        thicknesses.insert(0, self.ocean_model["depth"])
        dp = 0
        color_maps = ["Blues"] + ["Greys"]*(len(layers)-1)
        lev_values = [0.3] + np.linspace(.2,.8,len(layers)-1).tolist()
        for i in range(len(layers)):
            r = r_max - np.linspace(dp, layers[i]*d, N)
            t = np.linspace(t_min, t_max, N)
            ax.pcolormesh(t, r, lev_values[i]*np.ones((N, N)), vmin=0, vmax=1., cmap=color_maps[i])
            y, x = np.mean(r), np.mean(t)
            th = thicknesses[i]/1e3#efmte()
            if resistivities[i] > 1: txt = r"$\rho_E,\delta T_E\sim$ %d $\Omega m$,%d km"%(resistivities[i], th)
            else: txt = r"$\rho_E,\delta T_E\sim$ %.1f $\Omega m$,%d km"%(resistivities[i], th)
            #txt = r"$\rho_E,\delta T_E\sim $"
            if i == 0: txt = txt.replace("rho_E", "rho_O").replace("T_E", "T_O")
            else: txt = txt.replace("rho_E", "rho_E^%d"%i).replace("T_E", "T_E^%d"%i)
            ax.text(x-tag_offset, y, txt, color="firebrick", ha="center", va="center", fontdict=fontdict)
            dp = layers[i]*d
        r = r_max - np.linspace(dp, r_max, N)
        t = np.linspace(t_min, t_max, N)
        ax.pcolormesh(t, r, np.ones((N, N)), vmin=0, vmax=1., cmap="Greys")
        return
