"""models.py: Module is used to implement various sea-earth models and calculate TFs"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


# Import required packages
import bezpy
import matplotlib.pyplot as plt
import numpy as np
from scipy import constants as C

class LayeredOcean(object):
    """
    This class is to emulate electric and magnetic fields
    in a 1D ocean. We provide the 1D Earth conductivity model and
    one layer of ocean parameters (depth and resistivity) to emulate 
    1D ocean.
    
    Methods:
    --------
    update_resistivity_list: Update list of layers and resistivities
    calcZfloor: Caclulate Z for seafloor for given 1D Earth model and frequencies
    calcTF: Calculate various transfer functions
    plotTF: Plot various transfer functions with frequencies
    plot_resistivity: Plot resistivity depth plot
    plot_impedance_seafloor: Plot seafloor impedance with frequency
    
    Parameters:
    -----------
    model_name: Earth model name
    site: 1D Earth model
    ocean_parameters: Ocean parameters (depth, resistivity)
    depths: List of layer thickness (including ocean)
    resistivities: Resistivities in each layer (including ocean)
    freqs: Analysis frequency list in Hz
    """
    
    def __init__(self, model_name="SU1", ocean_parameters={"depth": 5e3, "rho":0.25}):
        """
        Initialize all model parameters.
        """
        self.model_name = model_name
        self.site = bezpy.mt.get_1d_site(model_name)
        self.ocean_parameters = ocean_parameters
        self.update_resistivity_list()
        return
    
    def update_resistivity_list(self, ocean_parameters=None):
        """
        Update list of layers and resistivities
        """
        if ocean_parameters is not None: self.ocean_parameters = ocean_parameters
        self.resistivities = np.array([self.ocean_parameters["rho"]] + self.site.resistivities.tolist())
        self.depths = np.array([self.ocean_parameters["depth"]] + (self.site.depths + self.ocean_parameters["depth"]).tolist())
        return
    
    def calcZfloor(self, flim=(1e-3, 1e3)):
        """
        Caclulate Z for seafloor for given 1D Earth model and frequencies
        """
        self.freqs = np.linspace(flim[0], flim[1], int(flim[1]/flim[0])+1)
        self.Zd = self.site.calcZ(self.freqs)[1,:]
        return
    
    def calcTF(self, kinds=["Ed2Ho", "Hd2Ho"], ax=None):
        """
        Calculate various transfer functions
        """
        omega = 2*C.pi*self.freqs
        sigma_s = 1/self.ocean_parameters["rho"]
        k2 = 1.j*omega*C.mu_0*sigma_s
        k = np.sqrt(k2)
        kd = k*self.ocean_parameters["depth"]
        Z = 1.j*omega*C.mu_0/k
        r = self.Zd/Z
        TFs = {}
        for kind in kinds:
            if kind=="Ed2Ho": tf = self.Zd / (np.cosh(kd) + r * np.sinh(kd))
            if kind=="Hd2Ho": tf = 1 / (np.cosh(kd) + r * np.sinh(kd))
            TFs[kind] = tf
        if ax is not None: self.plotTF(self.freqs, TFs, ax)
        return TFs
    
    def plotTF(self, freqs, TFs, ax):
        """
        Plot resistivity depth plot
        """
        ax.text(0.99, 1.05, r"$\rho_s (\Omega-m)$: %.2f"%(self.ocean_parameters["rho"]), 
            ha="right", va="center", transform=ax.transAxes)
        ax.text(1.05, 0.99, r"$D_{Ocean} (km)$: %d"%(self.ocean_parameters["depth"]/1e3), 
                ha="center", va="top", transform=ax.transAxes, rotation=90)
        ax.loglog(freqs, np.absolute(TFs["Ed2Ho"])*1e3, "r", lw=0.8, label=r"$\frac{E_d}{H_0}$")
        ax.loglog(freqs, np.absolute(TFs["Hd2Ho"]), "b", lw=0.8, label=r"$\frac{H_d}{H_0}$")
        ax.set_xlabel(r"$f_0$, (Hz)")
        ax.set_ylabel("Amplitude Ratio")
        ax.set_ylim(1e-10,1)
        ax.set_xlim(self.freqs[0],self.freqs[-1])
        ax.legend(loc=1)
        return
    
    def plot_resistivity(self, ax):
        """
        Plot resistivity depth plot
        """
        ax.text(0.99, 1.05, "Earth Model: %s"%self.site.name, ha="right", va="center", transform=ax.transAxes)
        ax.step(self.resistivities, np.insert(self.depths, 0, 0)/1000., color="b", lw=1.)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_ylabel("Depth (km)")
        ax.set_xlabel(r"Resistivity, $\rho$ ($\Omega-m$)")
        ax.set_ylim(1e-3, self.depths[-1]/1e3)
        ax.set_xlim(1e-3, 1e5)
        ax.invert_yaxis()
        return
    
    def plot_impedance_seafloor(self, ax):
        """
        Plot seafloor impedance with frequency
        """
        ax.text(0.99, 1.05, "Earth Model: %s"%self.site.name, ha="right", va="center", transform=ax.transAxes)
        mag, phase = np.abs(self.Zd), np.rad2deg(np.angle(self.Zd))
        ax.semilogx(self.freqs, mag, "r", lw=1.)
        ax.set_ylabel("|Z|=|a+jb|", fontdict={"color":"r"})
        ax.set_xlabel(r"$f_0$ (Hz)")
        ax = ax.twinx()
        ax.semilogx(self.freqs, phase, "b", lw=1.)
        ax.set_xlim(self.freqs[0], self.freqs[-1])
        ax.set_ylabel(r"$\theta(Z)=tan^{-1}(\frac{b}{a})$", fontdict={"color":"b"})
        return