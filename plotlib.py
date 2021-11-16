"""plotlib.py: Module is used to implement plotting functions"""

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

import datetime as dt

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

def plotZ(Z, freqs, ax, model_name=None, layer=None, 
          mag_ax=dict(color="r", label=r"$|Z|=|a+jb|$", lim=[1e-8, 1e0], lw=1., ls="-"), 
          pha_ax=dict(color="b", label=r"$\theta(Z)=\arctan(\frac{b}{a})$", lim=[0, 90], lw=1., ls="-")):
    """
    Plot impedance with frequency
    """
    model_name = "UD" if model_name is None else model_name
    layer = "UD" if layer is None else layer
    ax.text(0.99, 1.05, "Earth Model: %s (L=%s)"%(model_name, layer), ha="right", va="center", transform=ax.transAxes)
    mag, phase = np.abs(Z), np.rad2deg(np.angle(Z))
    ax.loglog(freqs, mag, mag_ax.color, lw=mag_ax.lw, ls=mag_ax.ls)
    ax.set_ylim(mag_ax.lim)
    ax.set_ylabel(mag_ax.label, fontdict={"color": mag_ax.color})
    ax.set_xlabel(r"$f_0$ (Hz)")
    ax = ax.twinx()
    ax.semilogx(self.freqs, phase, pha_ax.color, lw=pha_ax.lw, ls=pha_ax.ls)
    ax.set_ylim(pha_ax.lim)
    ax.set_ylabel(pha_ax.label, fontdict={"color": pha_ax.color})
    ax.set_xlim(self.freqs[0], self.freqs[-1])
    return ax

def plotTF(TFs, freqs, rs, th, ylims=[1e-2, 1e0]):
    """
    Plot transfer function frequency plot
    """
    if freqs is None: freqs = np.copy(self.freqs)
    ax.text(0.99, 1.05, r"$\rho_s (\Omega-m)$: %.2f"%(rs), ha="right", va="center", transform=ax.transAxes)
    ax.text(1.05, 0.99, r"$D_{Ocean} (km)$: %d"%(th/1e3), ha="center", va="top", transform=ax.transAxes, rotation=90)
    ax.loglog(freqs, np.absolute(TFs["Ef2Bs"]), "r", lw=0.8, label=r"$\left|\frac{E_f}{B_s}\right|$")
    ax.loglog(freqs, np.absolute(TFs["Ef2Bs"]), "b", lw=0.8, label=r"$\left|\frac{B_f}{B_s}\right|$")
    ax.set_xlabel(r"$f_0$, (Hz)")
    ax.set_ylabel("Amplitude Ratio")
    ax.set_ylim(ylims)
    ax.set_xlim(freqs[0],freqs[-1])
    ax.legend(loc=3)
    return

def create_pane(nrows=1, ncols=1, dpi=150, figsize=(3,3), wspace=0.2, hspace=0.2):
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, dpi=dpi, figsize=figsize)
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    return fig, axes


def plot_xy_magnetic_field(stns, frames, dpi=150, wspace=0.3, hspace=0.3):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    fig, axes = plt.subplots(nrows=len(stns), ncols=2, dpi=dpi, figsize=(12, 3*len(stns)))
    for i, stn in enumerate(stns):
        frame = frames[stn]
        ax = axes[i,0]
        #ax.tick_params(axis="y", which="both", colors="b")
        ax.set_xlim(frame.index.tolist()[0], frame.index.tolist()[-1]+dt.timedelta(minutes=1))
        ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
        ax.plot(frame.index, frame.X, "k", ls="-", lw=1.)
        ax.text(0.05, 0.9, stn.upper(), ha="left", va="center", transform=ax.transAxes, fontdict={"fontweight": "bold"})
        ax.set_ylabel(r"$B_x$, nT ($\times 10^{-9}$T)", fontdict={"color": "k"})
        if i == len(stns)-1: ax.set_xlabel("Time, UT")
        ax = axes[i,1]
        #ax.tick_params(axis="y", which="both", colors="r")
        ax.set_xlim(frame.index.tolist()[0], frame.index.tolist()[-1]+dt.timedelta(minutes=1))
        ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
        ax.plot(frame.index, frame.Y, "k", ls="-", lw=1.)
        ax.set_ylabel(r"$B_y$, nT ($\times 10^{-9}$T)", fontdict={"color": "k"})
        #ax.spines["left"].set_color("b")
        #ax.spines["right"].set_color("r")
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig("docs/Bxy.Field.png", bbox_inches="tight")
    return