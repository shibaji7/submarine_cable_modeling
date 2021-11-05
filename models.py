"""models.py: Module is used to implement various sea-earth models and calculate TFs"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
os.system("rm -rf `find -type d -name .ipynb_checkpoints`:")

# Import required packages
from bezpy.mt import Site1d
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import numpy as np
from scipy import constants as C
import pandas as pd
import cartopy.crs as ccrs

from multiprocessing import Pool
from fastkml import kml
from netCDF4 import Dataset

from math import radians, degrees, sin, cos, asin, acos, sqrt

def great_circle(lon1, lat1, lon2, lat2, R = 6371.):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    return R * ( acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)) )

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


    
def Ed2Ho(Z, Zd, kd):
    return Zd/(np.cosh(kd) + (Zd*np.sinh(kd)/Z))

def Hd2Ho(Z, Zd, kd):
    return 1./(np.cosh(kd) + (Zd*np.sinh(kd)/Z))
        
class LayeredOcean(object):
    """
    This class is to emulate electric and magnetic fields
    in a 1D ocean. We provide the 1D Earth conductivity model and
    one layer of ocean parameters (depth and resistivity) to emulate 
    1D ocean.
    
    Methods:
    --------
    update_resistivity_list: Update list of layers and resistivities
    calcTF: Calculate various transfer functions
    plotTF: Plot various transfer functions with frequencies
    plot_resistivity: Plot resistivity depth plot
    plot_impedance: Plot seafloor impedance with frequency
    
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
                 layers=None, flim=[1e-6, 1e0]):
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
        return
    
    def ini_tfs(self):
        """
        Create transfer (lambda) functions
        """
        self.functions = {"Ed2Ho": Ed2Ho, "Hd2Ho": Hd2Ho}
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
    
    def plot_resistivity(self, ax, params={"ylim":[0, 1e3], "xlim":[1e-3, 1e5], 
                                           "xscale":"log", "ylabel":"Depth (km)",
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
    
    def calcZ(self, layer=0, ocean=True):
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
            self.kd = k[0]*self.thicknesses[0]
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
            
            ## ###########################################
            ## Update 0th layer's impedance
            ## ###########################################
            if ocean:
                omega = 2*C.pi*self.freqs
                sigma_s = 1/self.resistivities[0]
                k2 = 1.j*omega*C.mu_0*sigma_s
                k = np.sqrt(k2)
                Z[0, :] = 1.j*omega*C.mu_0/k
            
            if freqs[0] == 0.: Z[:, 0] = 0.
            self.Z = np.copy(Z)
        
        else: nfreq = len(self.freqs)
        
        Z_output = np.zeros(shape=(4, nfreq), dtype=np.complex)
        Z_output[1, :] = self.Z[layer, :]
        Z_output[2, :] = -Z_output[1, :]
        return Z_output
    
    def calcTFRec(self, freq):
        """
        """
        resistivities = self.resistivities
        thicknesseses = self.depths
        omega = 2*np.pi*freq
        n = len(thicknesseses)
        
        omega = 2*np.pi*freq
        complex_factor = 1j*omega*C.mu_0
        k = np.sqrt(1j*omega*C.mu_0/resistivities)
        Z = np.zeros(shape=(n), dtype=np.complex)
        
        with np.errstate(divide="ignore", invalid="ignore"):
            Z[-1] = complex_factor/k[-1]
            
            r = np.zeros(shape=(n), dtype=np.complex)
            for i in range(n-2, -1, -1):
                r[i] = ((1-k[i]*Z[i]/complex_factor) /
                           (1+k[i]*Z[i+1]/complex_factor))
                Z[i] = (complex_factor*(1-r[i]*np.exp(-2*k[i]*thicknesseses[i])) /
                           (k[i]*(1+r[i]*np.exp(-2*k[i]*thicknesseses[i]))))
        o = self.functions["Ed2Ho"](Z[0], Z[1], k[0]*thicknesseses[0])
        return o
    
    def calcTF(self, kinds=["Ed2Ho", "Hd2Ho"], ax=None, ylims=[1e-2, 1e0], th=None):
        """
        Calculate various transfer functions
        """
        omega = 2*C.pi*self.freqs
        Zd = self.calcZ(1)[1, :]
        Z = self.calcZ(0)[1, :]
        TFs = {}
        
        omega = 2*C.pi*self.freqs
        sigma_s = 1/self.resistivities[0]
        k2 = 1.j*omega*C.mu_0*sigma_s
        k = np.sqrt(k2)
        kd = k*self.thicknesses[0] if th is None else k*th
        
        for kind in kinds:
            TFs[kind] = self.functions[kind](Z, Zd, kd)
        TFs["Ed2Bo"] = TFs["Ed2Ho"]*C.mu_0
        if ax is not None: self.plotTF(ax, TFs, ylims=ylims)
        return TFs
    
    def plotTF(self, ax, TFs, freqs=None, ylims=[1e-2, 1e0]):
        """
        Plot transfer function frequency plot
        """
        if freqs is None: freqs = np.copy(self.freqs)
        ax.text(0.99, 1.05, r"$\rho_s (\Omega-m)$: %.2f"%(self.resistivities[0]), 
            ha="right", va="center", transform=ax.transAxes)
        ax.text(1.05, 0.99, r"$D_{Ocean} (km)$: %d"%(self.thicknesses[0]/1e3), 
                ha="center", va="top", transform=ax.transAxes, rotation=90)
        ax.loglog(freqs, np.absolute(TFs["Ed2Ho"])*1e3, "r", lw=0.8, label=r"$\left|\frac{E_d}{B_o}\right|$")
        #ax.loglog(freqs, np.absolute(TFs["Ed2Bo"])*1e3, "r", lw=1.2, ls="--", label=r"$\frac{E_d}{B_o}$")
        ax.loglog(freqs, np.absolute(TFs["Hd2Ho"]), "b", lw=0.8, label=r"$\left|\frac{B_d}{B_o}\right|$")
        ax.set_xlabel(r"$f_0$, (Hz)")
        ax.set_ylabel("Amplitude Ratio")
        ax.set_ylim(ylims)
        ax.set_xlim(freqs[0],freqs[-1])
        ax.legend(loc=3)
        return
    
    def plotTFMagPhase(self, ax, freqs=None, ylims=[1e-2, 1e0], th=None):
        """
        Plot transfer function frequency plot
        """
        TFs = self.calcTF(th=th)
        if freqs is None: freqs = np.copy(self.freqs)
        ax.text(0.99, 1.05, r"$\rho_s (\Omega-m)$: %.2f"%(self.resistivities[0]), 
            ha="right", va="center", transform=ax.transAxes)
        ax.loglog(freqs, np.absolute(TFs["Ed2Ho"])*1e3, "r", lw=0.8)
        ax.set_xlabel(r"$f_0$, (Hz)")
        ax.set_ylabel(r"$\left|\frac{E_d}{B_o}\right|$", fontdict={"color":"r"})
        ax.set_ylim(ylims)
        ax.set_xlim(freqs[0],freqs[-1])
        ax = ax.twinx()
        ax.semilogx(freqs, 180*np.angle(TFs["Ed2Ho"])/np.pi, "b", lw=0.8)
        ax.set_ylim(0,90)
        ax.set_ylabel(r"$\theta(\frac{E_d}{B_o})$", fontdict={"color":"b"})
        return
    
    def plot_impedance(self, ax, layer=0):
        """
        Plot impedance with frequency for a leyer
        """
        ax.text(0.99, 1.05, "Earth Model: %s (L=%d)"%(self.model_name, layer), ha="right", va="center", transform=ax.transAxes)
        Z = self.calcZ(layer)[1, :]
        mag, phase = np.abs(Z), np.rad2deg(np.angle(Z))
        ax.loglog(self.freqs, mag, "r", lw=1.)
        ax.set_ylim(1e-8, 1e0)
        ax.set_ylabel(r"$|Z|=|a+jb|$", fontdict={"color":"r"})
        ax.set_xlabel(r"$f_0$ (Hz)")
        ax = ax.twinx()
        ax.semilogx(self.freqs, phase, "b", lw=1.)
        ax.set_ylim(0, 80)
        ax.set_xlim(self.freqs[0], self.freqs[-1])
        ax.set_ylabel(r"$\theta(Z)=\arctan(\frac{b}{a})$", fontdict={"color":"b"})
        return
    
class OceanFloor(object):
    """
    Ocean floor model based on parametric or emperical definations
    """
    
    def __init__(self, length_segment, depths, verbose=False):
        self.length_segment = length_segment
        self.depths = depths
        self.verbose = verbose
        return
    
    def plot_ocean_basin(self, ax):
        """
        Plot ocean basin based on segment and depths
        """
        ax.plot(self.length_segment, self.depths, "r")
        ax.axvline(0, ls="--", lw=1., color="k")
        ax.set_ylim(min(self.depths), max(self.depths))
        ax.set_xlim(min(self.length_segment), max(self.length_segment))
        ax.set_xlabel("Ocean Basin Segment/Length, km")
        ax.set_ylabel("Ocean Basin Depths, km")
        return
    
    def _run_rkck_(self, ls, d, model_name, flim):
        if self.verbose: print(" Calculate transfer functions for location, depths: %.2f, %.2f"%(ls, d))
        o = LayeredOcean(model_name=model_name, flim=flim)
        o.thicknesses[0] = d*1.e3
        o.calcZ()
        return o
    
    def create_basin_tf(self, model_name="BM1", flim=[1e-6, 1e0], nprocs=24):
        """
        Create a list of LayeredOcean for each segments
        """        
        if self.verbose: print(" Total processes (nprocs): %d(%d)"%(len(self.length_segment), nprocs))
        self.flim = flim
        inputs = [(ls, d, model_name, flim) for ls, d in zip(self.length_segment, np.abs(self.depths))]
        pool = Pool(nprocs)
        self.oceans = pool.starmap(self._run_rkck_, inputs)
        pool.close()
        pool.join()
        return
    
    def contour_plots(self, ax, param="Zd:mag", cmap="Greens", 
                      norm=mcolor.LogNorm(vmin=1e-7, vmax=1e-1), 
                      label=r"$|Z_d|$"):
        """
        Contour plot for a specific parameter
        """
        freqs, ls = np.linspace(self.flim[0], self.flim[1], int(self.flim[1]/self.flim[0])+1), self.length_segment
        pv = np.zeros((len(self.oceans), len(freqs)))
        for i, o in enumerate(self.oceans):
            if param=="Zd:mag": pv[i, :] = np.absolute(o.calcZ(1)[1, :])
            if param=="Zd:pha": pv[i, :] = np.angle(o.calcZ(1)[1, :])
            if param=="TF.Ed2Ho:mag": pv[i, :] = np.absolute(o.calcTF()["Ed2Ho"])
            if param=="TF.Ed2Ho:pha": pv[i, :] = np.angle(o.calcTF()["Ed2Ho"])
            if param=="TF.Hd2Ho:mag": pv[i, :] = np.absolute(o.calcTF()["Hd2Ho"])
            if param=="TF.Hd2Ho:pha": pv[i, :] = np.angle(o.calcTF()["Hd2Ho"])
        ax.contour(ls, freqs, pv.T, norm=norm, cmap=cmap)
        im = ax.contourf(ls, freqs, pv.T, cmap=cmap, norm=norm, vmax=1e-1, vmin=1e-7)
        cb = plt.gcf().colorbar(im, ax=ax, extend="max", shrink=0.5)
        cb.set_label(label)
        ax.set_yscale("log")
        ax.set_xlim(min(self.length_segment), max(self.length_segment))
        ax.set_xlabel("Ocean Basin Segment/Length, (km)")
        ax.set_ylim(freqs[0],freqs[-1])
        ax.set_ylabel(r"$f_0$, (Hz)")
        return
    
class SimpleParabolicOcean(OceanFloor):
    """
    Ocean floor model based on parabolic definations
    """
    
    def __init__(self, lseg=[-5,5], dmax=5., dseg=.2, verbose=False):
        self.lseg = lseg
        self.dmax = dmax
        self.npts = int((lseg[-1]-lseg[0])/dseg)+1
        x = np.linspace(lseg[0],lseg[-1],self.npts)
        a = np.max(x**2)/(dmax*4)
        y = x**2/(4*a)-dmax
        super().__init__(x, y, verbose=verbose)
        return

class CableDetails(object):
    """
    This class is dedicated to extract cable detailes form kml files
    """
    
    def __init__(self, cable_file="data/cablemap.info.kml", cables=["TAT-14"]):
        self.cable_file = cable_file
        self.cables = cables
        if ".kml" in cable_file:
            with open(cable_file, "rt") as f: 
                d = f.read()
                d = "\n".join(d.split("\n")[1:])
            self.root = kml.KML()
            self.root.from_string(d)
            self._extract_cable_details_()
        elif ".cable" in cable_file:
            self.cable_props = {}
            details = pd.read_csv(cable_file)
            for c in self.cables:
                self.cable_props[c] = {"name": c, "lats": [], "lons": [], "alts": [], "A->E": {}}
                self.cable_props[c]["lons"].extend(details.lons)
                self.cable_props[c]["lats"].extend(details.lats)
                self.cable_props[c]["alts"].extend([0]*len(details))
                self.cable_props[c]["A->E"]["lats"] = details.lats.tolist()
                self.cable_props[c]["A->E"]["lons"] = details.lons.tolist()
                self.cable_props[c]["A->E"]["thickness"] = (1e-3*details.thickness).tolist()
        else: print(f"System is not able to work with file type .{cable_file.split('.')[-1]}")
        return
    
    def _extract_cable_details_(self):
        folders = list(list(list(self.root.features())[0].features())[1].features())
        self.cable_props = {}
        for f in folders:
            if f.name in self.cables: 
                cfeatures = list(f.features())
                self.cable_props[f.name] = {"name": f.name, "lats": [], "lons": [], "alts": []}
                coords = list(cfeatures[0].geometry.coords)
                for c in coords:
                    self.cable_props[f.name]["lons"].append(c[0])
                    self.cable_props[f.name]["lats"].append(c[1])
                    self.cable_props[f.name]["alts"].append(c[2])
        return
    
    def get_ocean_depth_data(self, fname="data/LITHO1.0.nc"):
        nc = Dataset(fname)
        latitude, longitude = nc.variables["latitude"][:], nc.variables["longitude"][:]
        ocean_thickness = nc.variables["water_bottom_depth"][:] - nc.variables["water_top_depth"][:]
        thickness = []
        for c in self.cables:
            cprops = self.cable_props[c]
            clats, clons = np.array(cprops["lats"]), np.array(cprops["lons"])
            for lat, lon in zip(clats, clons):
                i, j = np.argmin(np.abs(latitude-lat)), np.argmin(np.abs(longitude-lon))
                thickness.append(ocean_thickness[i, j])
            first_index, last_index = None, None
            for _i, t in enumerate(thickness):
                if (not np.isnan(t)) and (first_index is None): first_index = _i
                if (np.isnan(t)) and (first_index is not None) and (last_index is None): last_index = _i
            self.cable_props[c]["A->E"] = {}
            self.cable_props[c]["A->E"]["lats"] = clats[first_index-1:last_index+1]
            self.cable_props[c]["A->E"]["lons"] = clons[first_index-1:last_index+1]
            self.cable_props[c]["A->E"]["thickness"] = thickness[first_index-1:last_index+1]
        nc.close()
        return
    
    def create_basin_tf(self, model_name="BM1", flim=[1e-6, 1e0]):
        for c in self.cables:
            self.cable_props[c]["A->E"]["oceans"] = []
            x = [(lat, lon, th) for lat, lon, th in zip(self.cable_props[c]["A->E"]["lats"], self.cable_props[c]["A->E"]["lons"],
                                                       self.cable_props[c]["A->E"]["thickness"])]
            for d in self.cable_props[c]["A->E"]["thickness"]:
                o = LayeredOcean(model_name=model_name, flim=flim)
                o.thicknesses[0] = d*1.e3
                o.calcZ()
                self.cable_props[c]["A->E"]["oceans"].append(o)
        return
    
    def plot_cable_footprint(self, loc=None, ax=None):
        if ax is None: ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        for c in self.cables:
            ax.plot(self.cable_props[c]["A->E"]["lons"],self.cable_props[c]["A->E"]["lats"], linewidth=1.,
                   transform=ccrs.PlateCarree(), color="b")
        ax.set_extent([-100, 30, 0, 80], crs=ccrs.PlateCarree())
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              linewidth=0.5, color="gray", alpha=0.3, linestyle="--")
        if loc is not None: ax.scatter([loc[0]], [loc[1]], s=5, marker="D", color="r", transform=ccrs.PlateCarree())
        return
    
    def plot_TF_location(self, ix=3, c="TAT-8"):
        fig = plt.figure(dpi=150, figsize=(4,8))
        ax = fig.add_subplot(211, projection=ccrs.PlateCarree())
        loc = (self.cable_props[c]["A->E"]["lons"][ix], self.cable_props[c]["A->E"]["lats"][ix])
        self.plot_cable_footprint(loc, ax)
        o = self.cable_props[c]["A->E"]["oceans"][ix]
        ax = fig.add_subplot(212)
        _ = o.calcTF(ax=ax)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        return
    
    
    
