import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import numpy as np
import pandas as pd
from scipy import constants as C

from scubas.datasets import Site
from scubas.utils import fft, ifft

def calcZ(site, freqs):
    omega = 2 * C.pi * freqs
    sigma_s = 1 / site.layers[0].resistivity
    k2 = 1.0j * omega * C.mu_0 * sigma_s
    k = np.sqrt(k2)
    Zo = (1.0j * omega * C.mu_0 / k) / (C.mu_0 / 1.0e-3)
    return Zo

def calcTF(site, freqs):
    """
    Calculate the transfer functions.
    """
    Zo = calcZ(site, freqs)
    Zd = Zo
    omega = 2 * C.pi * freqs
    sigma_s = 1 / site.layers[0].resistivity
    k2 = 1.0j * omega * C.mu_0 * sigma_s
    k = np.sqrt(k2)
    kd = k * 0

    func = Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Zo))
    tf = pd.DataFrame()
    tf["freq"], tf["E2B"] = freqs, func
    return tf

def calcTFx(site, freqs):
    """
    Calculate the transfer functions.
    """
    Zo = calcZ(site, freqs)
    Zd = site.calcZ(freqs)[1]
    omega = 2 * C.pi * freqs
    sigma_s = 1 / site.layers[0].resistivity
    k2 = 1.0j * omega * C.mu_0 * sigma_s
    k = np.sqrt(k2)
    kd = k * 0

    func = Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Zo))
    tf = pd.DataFrame()
    tf["freq"], tf["E2B"] = freqs, func
    return tf

def calcAppResist(site, freqs):
    tf = calcTFx(site, freqs)
    tf["omega"] = 2*np.pi*freqs
    tf["rho"] = np.abs(tf.E2B)**2/(2*np.pi*freqs*C.mu_0)
    tf["sqrt_T"] = np.sqrt(1/freqs)
    tf["skin_depth"] = 1/np.abs(np.sqrt(1j*tf["omega"]*C.mu_0*site.layers[0].conductivity))
    tf["nskin_depth"] = site.layers[0].thickness/tf["skin_depth"]
    return tf

sites = [
    Site.init(
        conductivities=[1/3, 1/3000],
        thicknesses=[10000, np.inf],
        names=[
            "Sediments",
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper Case A",
        site_name="Two layer Earth model",
    ),
    Site.init(
        conductivities=[1/3000, 1/3],
        thicknesses=[100000, np.inf],
        names=[
            "Sediments",
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper Case B",
        site_name="Two layer Earth model",
    )
]

flim, M = [1e-10, 1e0], 500
freqs = np.linspace(flim[0], flim[1], int(flim[1] / (M*flim[0])) + 1)
tfs = [
    calcAppResist(sites[0], freqs),
    calcAppResist(sites[1], freqs),
]

fig = plt.figure(dpi=300, figsize=(3, 3))
ax = fig.add_subplot(111)
ax.loglog(tfs[0].skin_depth/1e3, np.abs(tfs[0].rho), "r", lw=1.0, ls="-", label="Case A")
ax.loglog(tfs[1].skin_depth/1e3, np.abs(tfs[1].rho), "b", lw=1.0, ls="-", label="Case B")
ax.legend(loc=1)
ax.set_ylabel(r"$\rho_{ef}=\frac{|Z_{ef}|^2}{\omega\mu_0}$ [$\Omega-m$]")
ax.set_xlabel(r"$p_1=|\frac{1}{\sqrt{i\omega\mu_0\sigma}}|$ [$km$]")
ax.set_xlim(1e0, 1e4)
ax.set_ylim(1e6, 1e10)
fig.savefig("rho.png", bbox_inches="tight")

fig = plt.figure(dpi=300, figsize=(3, 3))
ax = fig.add_subplot(111)
ax.loglog(tfs[0].sqrt_T, np.abs(tfs[0].rho), "r", lw=1.0, ls="-", label="Case A")
ax.loglog(tfs[1].sqrt_T, np.abs(tfs[1].rho), "b", lw=1.0, ls="-", label="Case B")
ax.legend(loc=1)
ax.set_ylabel(r"$\rho_{ef}=\frac{|Z_{ef}|^2}{\omega\mu_0}$ [$\Omega-m$]")
ax.invert_xaxis()
ax.set_xlabel(r"$\sqrt{T}$ [$s^{1/2}$]")
ax.set_xlim([1e0, 1e4])
ax.set_ylim(1e6, 1e10)
fig.savefig("rho1.png", bbox_inches="tight")

fig = plt.figure(dpi=300, figsize=(3, 3))
ax = fig.add_subplot(111)
ax.loglog(tfs[0].nskin_depth, np.abs(tfs[0].rho), "r", lw=1.0, ls="-", label="Case A")
ax.loglog(tfs[1].nskin_depth, np.abs(tfs[1].rho), "b", lw=1.0, ls="-", label="Case B")
ax.legend(loc=1)
ax.set_ylabel(r"$\rho_{ef}=\frac{|Z_{ef}|^2}{\omega\mu_0}$ [$\Omega-m$]")
ax.set_xlabel(r"$\frac{\tau_1}{p_1}$")
ax.set_xlim(1e-4, 1e1)
ax.set_ylim(1e6, 1e10)
fig.savefig("rho3.png", bbox_inches="tight")
