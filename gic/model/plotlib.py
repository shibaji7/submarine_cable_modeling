"""plotlib.py: Module is used to plotting tools"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import matplotlib as mpl
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.style.use(["science", "ieee"])

# import sys

# sys.path.extend(["py/", "py/config/"])
from .utils import *
import numpy as np
from scipy.stats import pearsonr


class Summary(object):
    """
    Summary plots for the analysis data
    """

    def __init__(self, nrows=1, ncols=1, dpi=180, size=(5, 5)):
        self.nrows = nrows
        self.ncols = ncols
        self.dpi = dpi
        self.size = (size[0] * self.ncols, size[1] * self.nrows)
        self.fig = plt.figure(dpi=dpi, figsize=size)
        self.fnum = 0
        return

    def add_axes(self):
        self.fnum += 1
        ax = self.fig.add_subplot(self.nrows, self.ncols, self.fnum)
        return ax

    def save(self, fname):
        self.fig.subplots_adjust(wspace=0.7, hspace=0.7)
        self.fig.savefig(fname, bbox_inches="tight")
        return

    def close(self):
        plt.close()
        return


class BfieldSummary(Summary):
    """
    B-Field summary plot
    """

    def __init__(self, nrows=2, ncols=2, dpi=180, size=(5, 5)):
        super().__init__(nrows, ncols, dpi, size)
        return

    def add_Bfield_Seq(self, B, E):
        """
        Add synthetic B-field data
        """
        ylim = [(int(np.min(B.X / 10)) - 1) * 10, (int(np.max(B.X / 10)) + 1) * 10]
        xlim = [np.min(B.dTime / 3600.0), np.max(B.dTime / 3600.0)]
        ax = self.add_axes()
        ax.plot(B.dTime / 3600.0, B.X, ls="-", lw=0.8)
        ax.set_xlabel("Time, Hours")
        ax.set_ylabel("B-Field, nT")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        ax = ax.twinx()
        ax.plot(E.dTime / 3600.0, E.X, color="r", ls="-", lw=0.8)
        ax.set_ylabel("E-Field, mv/km", color="r")
        ylim = [(int(np.min(E.X / 10)) - 1) * 10, (int(np.max(E.X / 10)) + 1) * 10]
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        return ax

    def add_Es(self, Ea, En):
        ax = self.add_axes()
        ax.plot(Ea, En, "ko", ms=0.1, alpha=0.4)
        ax.set_xlabel(r"$E^{anl}(t)$")
        ax.set_ylabel(r"$E^{fft}(t)$")
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="r", ls="--", lw=0.8)
        r, _ = pearsonr(Ea, En)
        ax.text(
            0.1,
            0.9,
            r"$\rho=$%.10f" % r,
            va="center",
            ha="left",
            transform=ax.transAxes,
        )
        # ax.set_xlim([-15, 15])
        # ax.set_ylim([-15, 15])
        return ax

    def add_TK_param(self, tf):
        ax = self.add_axes()
        ax.semilogx(tf.freq, np.abs(tf.E2B), "k", lw=0.4)
        ax.set_xlim(1e-4, 1e-2)
        ax.set_ylabel(r"$|X(f)|$")
        ax.set_xlabel(r"$f_0$, Hz")
        ax = ax.twinx()
        ax.semilogx(tf.freq, np.angle(tf.E2B, deg=True), "r", lw=0.4)
        ax.set_xlim(1e-4, 1e-2)
        ax.set_ylabel(r"$\theta[X(f)]$", color="r")
        return ax


class AnalysisSummary(Summary):
    """
    Simulation summary plots
    """

    def __init__(self, nrows=2, ncols=2, dpi=180, size=(5, 5)):
        super().__init__(nrows, ncols, dpi, size)
        return


def potential_along_section(
    V, x, pname, sec=None, Vi=None, Vk=None, Z=None, Y=None, gma=None, Z0=None
):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize": 12, "font.size": 12})
    fig, axes = plt.subplots(
        nrows=1, ncols=1, dpi=150, figsize=(6, 3), sharex="all", sharey="all"
    )
    ax = axes
    ax.set_ylabel("Voltage, V")
    ax.set_xlabel("Cable Length, km")
    ax.plot(x, V, "k", lw=0.8, ls="-")
    if Z is not None:
        Z *= 1e3
    if Y is not None:
        Y *= 1e3
    if gma is not None:
        gma *= 1e3
    txt = ""
    if sec is not None:
        txt += "Along: Bin%02d\n" % (sec)
    if (Vi is not None) and (Vk is not None):
        txt += r"$V_i,V_k\sim %.1f V, %.1f V$" % (Vi, Vk) + "\n"
    if (Z is not None) and (Y is not None):
        txt += (
            r"$Z,Y\sim$ %s $\Omega/km$, %s $\mho/km$" % (frexp102str(Z), frexp102str(Y))
            + "\n"
        )
    if (gma is not None) and (Z0 is not None):
        txt += (
            r"$\gamma,Z_0\sim$ %s /km, %s $\Omega$"
            % (frexp102str(gma), frexp102str(Z0))
            + "\n"
        )
    txt += "L=%d km" % np.max(x)
    ax.text(
        0.05, 0.95, txt, ha="left", va="top", transform=ax.transAxes, fontsize="small"
    )
    ax.set_xlim(x[0], x[-1])
    fig.savefig(pname, bbox_inches="tight")
    return


def cable_potential(V, x, pname):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize": 12, "font.size": 12})
    fig, axes = plt.subplots(
        nrows=1, ncols=1, dpi=150, figsize=(6, 3), sharex="all", sharey="all"
    )
    ax = axes
    ax.set_ylabel("Voltage, V")
    ax.set_xlabel("Cable Length, km")
    ax.plot(x, V, "k", lw=0.8, ls="-")
    ax.set_xlim(x[0], x[-1])
    txt = ""
    ax.text(
        0.05, 0.95, txt, ha="left", va="top", transform=ax.transAxes, fontsize="small"
    )
    fig.savefig(pname, bbox_inches="tight")
    return
