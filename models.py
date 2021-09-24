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
from bezpy.mt import Site1d
import matplotlib.pyplot as plt
import numpy as np
from scipy import constants as C

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


def read_1d_usgs_profile(fname):
    """Reads in a USGS conductivity profile.
    These are downloaded from:
    https://geomag.usgs.gov/conductivity/index.php
    Note that they are of a specific format, with thicknesses and conductivity
    listed. This should be adaptable to other 1d profiles from different locations
    with minor modifications to return a new 1d site object.
    """

    conductivities = []
    thicknesses = []
    # Just keep the last part of the file name
    profile_name = fname.strip(".txt").split("_")[-1]

    with open(fname, "r") as f:
        for line in f:
            if line[0] == "*":
                continue
            # Moved past all comment lines
            # First line is supposed to be the number of layers
            num_layers = int(line.split()[0])
            f.readline()  # Spaces between each set of points

            for _ in range(num_layers):
                # Alternates conductivity/depth
                conductivities.append(float(f.readline().split()[0]))
                thicknesses.append(float(f.readline().split()[0]))
                f.readline()
            conductivities.append(float(f.readline().split()[0]))

            # Done with the file, so create the site and return it
            site = Site1d(name=profile_name,
                          thicknesses=thicknesses,
                          resistivities=[1./x for x in conductivities])
            return site

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
    depths: List of layer thicknesses (including ocean)
    resistivities: Resistivities in each layer (including ocean)
    freqs: Analysis frequency list in Hz
    """
    
    def __init__(self, model_name="BM1", ocean_parameters=None, 
                 layers=None, flim=[1e-3, 1e3], layer_interaction=False):
        """
        Initialize all model parameters.
        """
        self.model_name = model_name
        self.site = read_1d_usgs_profile("data/ocean_model_%s.txt"%model_name)
        self.ocean_parameters = ocean_parameters
        self.layers = layers
        self.freqs = np.linspace(flim[0], flim[1], int(flim[1]/flim[0])+1)
        self.ini_tfs()
        self.update_resistivity_list()
        self.layer_interaction = layer_interaction
        return
    
    def ini_tfs(self):
        """
        Create transfer (lambda) functions
        """
        self.functions = {}
        self.functions["Ed2Ho"] = lambda Z, Zd, kd: Zd/(np.cosh(kd) + (Zd*np.sinh(kd)/Z))
        self.functions["Hd2Ho"] = lambda Z, Zd, kd: 1./(np.cosh(kd) + (Zd*np.sinh(kd)/Z))
        return
    
    @property
    def depths(self):
        """The depth of the profiles."""
        return np.cumsum(self.thicknesses)
    
    def update_resistivity_list(self):
        """
        Update list of layers and resistivities
        """
        if self.layers is not None:
            self.resistivities, self.thicknesses = self.layers["rho"], self.layers["thickness"]
        elif self.ocean_parameters is not None:
            self.resistivities = np.array([self.ocean_parameters["rho"]] + self.site.resistivities.tolist())
            self.thicknesses = np.array([self.ocean_parameters["thickness"]] + self.site.thicknesses.tolist())
        else: self.resistivities, self.thicknesses = self.site.resistivities, self.site.thicknesses
        return
    
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
        N = params["N"] if "N" in params.keys() else 1001
        t_min = params["t_min"] if "t_min" in params.keys() else 60.
        t_max = params["t_max"] if "t_max" in params.keys() else 120.
        r_min =  params["r_min"] if "r_min" in params.keys() else 0.
        tag_offset = params["tag_offset"] if "tag_offset" in params.keys() else 0.4
        thick = params["thick"] if "thick" in params.keys() else 2
        layers = np.arange(1, len(self.thicknesses)+1)
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
        
        conductance = 1./self.resistivities
        dp = 0
        color_maps = ["Blues"] + ["Greys"]*(len(layers)-1)
        lev_values = [0.3] + np.linspace(.2,.8,len(layers)-1).tolist()
        for i in range(len(layers)):
            r = r_max - np.linspace(dp, layers[i]*d, N)
            t = np.linspace(t_min, t_max, N)
            ax.pcolormesh(t, r, lev_values[i]*np.ones((N, N)), vmin=0, vmax=1., cmap=color_maps[i])
            y, x = np.mean(r), np.mean(t)
            th = efmte(self.thicknesses[i]/1e3)
            txt = r"$\sigma_e,\delta t_e\sim %.2f,%s$"%(conductance[i], th)
            if i == 0: txt = txt.replace("sigma_e", "sigma_o").replace("t_e", "t_o")
            ax.text(x-tag_offset, y, txt, color="firebrick", ha="center", va="center", fontdict=fontdict)
            dp = layers[i]*d
        r = r_max - np.linspace(dp, r_max, N)
        t = np.linspace(t_min, t_max, N)
        ax.pcolormesh(t, r, np.ones((N, N)), vmin=0, vmax=1., cmap="Greys")
        return
    
    def plot_resistivity(self, ax, params={"ylim":[1e-3, 1e3], "xlim":[1e-3, 1e5], 
                                           "yscale":"log", "xscale":"log", "ylabel":"Depth (km)",
                                           "xlabel":r"Resistivity, $\rho$ ($\Omega-m$)"}):
        """
        Plot resistivity depth plot
        """
        ax.text(0.99, 1.05, "Earth Model: %s"%self.model_name, ha="right", va="center", transform=ax.transAxes)
        ax.step(self.resistivities, np.insert(self.depths, 0, 0)/1000., color="k", lw=1.)
        ax.step([self.resistivities[0],self.resistivities[0]], np.insert(self.depths, 0, 0)[:2]/1000., color="b", lw=1.)
        if "yscale" in params.keys(): ax.set_yscale(params["yscale"])
        if "xscale" in params.keys(): ax.set_xscale(params["xscale"])
        if "ylabel" in params.keys(): ax.set_ylabel(params["ylabel"])
        if "xlabel" in params.keys(): ax.set_xlabel(params["xlabel"])
        if "ylim" in params.keys(): ax.set_ylim(params["ylim"][0], self.depths[-1]/1e3)
        if "xlim" in params.keys(): ax.set_xlim(params["xlim"][0], params["xlim"][1])
        ax.invert_yaxis()
        return
    
    def calcZfloor(self):
        """
        Caclulate Z for seafloor for given 1D Earth model and frequencies
        """
        self.Zf = self.calcZ(0)[1, :]
        return
    
    def calcZsurface(self):
        """
        Caclulate Z for sea surface for given 1D Earth model and frequencies
        """
        self.Zs = self.calcZ(1)[1, :]
        return
    
    def calcZ(self, layer=0):
        """
        Caclulate Zs for given 1D Earth model and frequencies
        """
        if not hasattr(self, "Z"):
            
            freqs = np.copy(self.freqs)
            resistivities = self.resistivities
            thicknesseses = self.depths
        
            n = len(resistivities)
            nfreq = len(freqs)

            omega = 2*np.pi*freqs
            complex_factor = 1j*omega*C.mu_0

            k = np.sqrt(1j*omega[np.newaxis, :]*C.mu_0/resistivities[:, np.newaxis])
            Z = np.zeros(shape=(n, nfreq), dtype=np.complex)
            # DC frequency produces divide by zero errors
            with np.errstate(divide="ignore", invalid="ignore"):
                Z[-1, :] = complex_factor/k[-1, :]

                r = np.zeros(shape=(n, nfreq), dtype=np.complex)
                for i in range(n-2, -1, -1):
                    r[i, :] = ((1-k[i, :]*Z[i+1, :]/complex_factor) /
                               (1+k[i, :]*Z[i+1, :]/complex_factor))
                    Z[i, :] = (complex_factor*(1-r[i, :]*np.exp(-2*k[i, :]*thicknesseses[i])) /
                               (k[i, :]*(1+r[i, :]*np.exp(-2*k[i, :]*thicknesseses[i]))))
            if freqs[0] == 0.: Z[:, 0] = 0.
            self.Z = np.copy(Z)
        else: nfreq = len(self.freqs)
        Z_output = np.zeros(shape=(4, nfreq), dtype=np.complex)
        Z_output[1, :] = self.Z[layer, :]*(1.e-3/C.mu_0)
        Z_output[2, :] = -Z_output[1, :]
        return Z_output
    
    def calcTF(self, kinds=["Ed2Ho", "Hd2Ho"], ax=None):
        """
        Calculate various transfer functions
        """
        omega = 2*C.pi*self.freqs
        Zd = self.calcZ(1)[1, :]
        if self.layer_interaction:
            Z, kd = self.calcZ(0)[1, :], (np.sqrt(1j*omega[np.newaxis, :]*C.mu_0/
                                                  self.resistivities[:, np.newaxis])[0])*self.thicknesses[0]
        else:
            sigma_s = 1/self.resistivities[0]
            k2 = 1.j*omega*C.mu_0*sigma_s
            k = np.sqrt(k2)
            kd = k*self.thicknesses[0]
            Z = 1.j*omega*C.mu_0/k
        TFs = {}
        for kind in kinds:
            TFs[kind] = self.functions[kind](Z, Zd, kd)
        if ax is not None: self.plotTF(ax, TFs)
        return TFs
    
    def plotTF(self, ax, TFs, freqs=None, ylims=[1e-10, 1]):
        """
        Plot transfer function frequency plot
        """
        if freqs is None: freqs = np.copy(self.freqs)
        ax.text(0.99, 1.05, r"$\rho_s (\Omega-m)$: %.2f"%(self.resistivities[0]), 
            ha="right", va="center", transform=ax.transAxes)
        ax.text(1.05, 0.99, r"$D_{Ocean} (km)$: %d"%(self.thicknesses[0]/1e3), 
                ha="center", va="top", transform=ax.transAxes, rotation=90)
        ax.loglog(freqs, np.absolute(TFs["Ed2Ho"])*1e3, "r", lw=0.8, label=r"$\frac{E_d}{H_0}$")
        ax.loglog(freqs, np.absolute(TFs["Hd2Ho"]), "b", lw=0.8, label=r"$\frac{H_d}{H_0}$")
        ax.set_xlabel(r"$f_0$, (Hz)")
        ax.set_ylabel("Amplitude Ratio")
        ax.set_ylim(ylims)
        ax.set_xlim(freqs[0],freqs[-1])
        ax.legend(loc=1)
        return
    
    def plot_impedance(self, ax, layer=0):
        """
        Plot impedance with frequency for a leyer
        """
        ax.text(0.99, 1.05, "Earth Model: %s (L=%d)"%(self.model_name, layer), ha="right", va="center", transform=ax.transAxes)
        if (not self.layer_interaction) and layer == 0: 
            omega = 2*C.pi*self.freqs
            sigma_s = 1/self.resistivities[0]
            k2 = 1.j*omega*C.mu_0*sigma_s
            k = np.sqrt(k2)
            Z = 1.j*omega*C.mu_0/k
        else: Z = self.calcZ(layer)[1, :]
        mag, phase = np.abs(Z), np.rad2deg(np.angle(Z))
        ax.semilogx(self.freqs, mag, "r", lw=1.)
        ax.set_ylabel("|Z|=|a+jb|", fontdict={"color":"r"})
        ax.set_xlabel(r"$f_0$ (Hz)")
        ax = ax.twinx()
        ax.semilogx(self.freqs, phase, "b", lw=1.)
        ax.set_xlim(self.freqs[0], self.freqs[-1])
        ax.set_ylabel(r"$\theta(Z)=tan^{-1}(\frac{b}{a})$", fontdict={"color":"b"})
        return