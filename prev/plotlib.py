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
import numpy as np

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

def plotTF(ax, TFs, freqs, rs, th, ylims=[1e-2, 1e0]):
    """
    Plot transfer function frequency plot
    """
    if freqs is None: freqs = np.copy(self.freqs)
    ax.text(0.99, 1.05, r"$\rho_s (\Omega-m)$: %.2f"%(rs), ha="right", va="center", transform=ax.transAxes)
    ax.text(1.05, 0.99, r"$D_{Ocean} (km)$: %.2f"%(th/1e3), ha="center", va="top", transform=ax.transAxes, rotation=90)
    ax.loglog(freqs, np.absolute(TFs["Ef2Bs"]), "r", lw=0.8, label=r"$\left|\frac{E_f}{B_s}\right|$")
    ax.loglog(freqs, np.absolute(TFs["Bf2Bs"]), "b", lw=0.8, label=r"$\left|\frac{B_f}{B_s}\right|$")
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

def plot_Bxy(stn, frame, fname, dpi=150, wspace=0.2, hspace=0.1):
    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=(5,3))
    ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
    ax.set_ylabel(r"$B_{x}$, nT", fontdict={"color": "b"})
    ax.set_xlabel(r"Time, UT", fontdict={"color": "k"})
    ax.text(0.05, 0.9, stn.upper(), ha="left", va="center", transform=ax.transAxes, 
            fontdict={"fontweight": "bold"})
    col = "darkblue"
    ax.plot(frame.index, frame.X, col, ls="-", lw=1.)
    ax.spines["left"].set_color(col)
    ax.tick_params(axis="y", which="both", colors=col)
    ax.yaxis.label.set_color(col)
    ax = ax.twinx()
    ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
    col = "darkred"
    ax.plot(frame.index, frame.Y, col, ls="-", lw=1.)
    ax.set_ylabel(r"$B_{y}$, nT", fontdict={"color": col})
    ax.spines["right"].set_color(col)
    ax.tick_params(axis="y", which="both", colors=col)
    ax.yaxis.label.set_color(col)
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig(fname, bbox_inches="tight")
    return

def plot_Exy(stn, frame, fname, dpi=150, wspace=0.2, hspace=0.1):
    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=(5,3))
    ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
    ax.set_ylabel(r"$E_{x}$, $mV.km^{-1}$", fontdict={"color": "b"})
    ax.set_xlabel(r"Time, UT", fontdict={"color": "k"})
    ax.text(0.05, 0.9, stn.upper(), ha="left", va="center", transform=ax.transAxes, 
            fontdict={"fontweight": "bold"})
    col = "darkblue"
    ax.plot(frame.index, frame.X, col, ls="-", lw=1.)
    ax.spines["left"].set_color(col)
    ax.tick_params(axis="y", which="both", colors=col)
    ax.yaxis.label.set_color(col)
    ax = ax.twinx()
    ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
    col = "darkred"
    ax.plot(frame.index, frame.Y, col, ls="-", lw=1.)
    ax.set_ylabel(r"$E_{y}$, $mV.km^{-1}$", fontdict={"color": col})
    ax.spines["right"].set_color(col)
    ax.tick_params(axis="y", which="both", colors=col)
    ax.yaxis.label.set_color(col)
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig(fname, bbox_inches="tight")
    return

def plot_BExy(stn, Bframe, Eframe, fname, dpi=150, wspace=0.2, hspace=0.1):
    fig, axes = plt.subplots(nrows=2, ncols=1, dpi=dpi, figsize=(4,5))
    ax = axes[0]
    ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
    ax.set_ylabel(r"$B_{x}$, nT", fontdict={"color": "k"})
    ax.text(0.05, 0.9, stn.upper(), ha="left", va="center", transform=ax.transAxes, 
            fontdict={"fontweight": "bold"})
    col = "darkblue"
    ax.plot(Bframe.index, Bframe.X, col, ls="-", lw=1.)
    ax.spines["left"].set_color(col)
    ax.tick_params(axis="y", which="both", colors=col)
    ax.yaxis.label.set_color(col)
    ax = ax.twinx()
    ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
    col = "darkred"
    ax.plot(Eframe.index, Eframe.X, col, ls="-", lw=1.)
    ax.set_ylabel(r"$E_{x}$, $mV.km^{-1}$", fontdict={"color": col})
    ax.spines["right"].set_color(col)
    ax.tick_params(axis="y", which="both", colors=col)
    ax.yaxis.label.set_color(col)
    
    ax = axes[1]
    ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
    ax.set_ylabel(r"$B_{y}$, nT", fontdict={"color": "k"})
    ax.set_xlabel(r"Time, UT", fontdict={"color": "k"})
    col = "darkblue"
    ax.plot(Bframe.index, Bframe.Y, col, ls="-", lw=1.)
    ax.spines["left"].set_color(col)
    ax.tick_params(axis="y", which="both", colors=col)
    ax.yaxis.label.set_color(col)
    ax = ax.twinx()
    ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
    col = "darkred"
    ax.plot(Eframe.index, Eframe.Y, col, ls="-", lw=1.)
    ax.set_ylabel(r"$E_{y}$, $mV.km^{-1}$", fontdict={"color": col})
    ax.spines["right"].set_color(col)
    ax.tick_params(axis="y", which="both", colors=col)
    ax.yaxis.label.set_color(col)
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig(fname, bbox_inches="tight")
    return

def plot_Exy_Stack(frames, fname, dpi=150, wspace=0.2, hspace=0.1):
    fig, axes = plt.subplots(nrows=len(frames), ncols=1, dpi=dpi, figsize=(4,2.5*len(frames)))
    for ax, k in zip(axes, frames.keys()):
        frame = frames[k]
        ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
        ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
        col = "darkblue"
        ax.plot(frame.index, frame.X, col, ls="-", lw=1.)
        ax.set_ylabel(r"$E_{x}$, $mV.km^{-1}$", fontdict={"color": col})
        ax.spines["left"].set_color(col)
        ax.tick_params(axis="y", which="both", colors=col)
        ax.yaxis.label.set_color(col)
        ax = ax.twinx()
        ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
        ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
        col = "darkred"
        ax.plot(frame.index, frame.Y, col, ls="-", lw=1.)
        ax.set_ylabel(r"$E_{y}$, $mV.km^{-1}$", fontdict={"color": col})
        ax.spines["right"].set_color(col)
        ax.tick_params(axis="y", which="both", colors=col)
        ax.yaxis.label.set_color(col)
        pass
    ax.set_xlabel(r"Time, UT", fontdict={"color": "k"})
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig(fname, bbox_inches="tight")
    return

def plot_xy_magnetic_field_oneplot(stns, frames, dpi=150, wspace=0.2, hspace=0.1):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    fig, axes = plt.subplots(nrows=1, ncols=2, dpi=dpi, figsize=(12, len(stns)))
    multiplier, colors = [1, 0, -1], ["r", "k", "b"]
    base=1000
    for i, stn in enumerate(stns):
        frame = frames[stn]
        ax = axes[0]
        ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
        ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
        ax.plot(frame.index, (base*multiplier[i])+frame.X-np.mean(frame.X), colors[i], ls="-", lw=1., label=stn.upper())
        ax.set_ylabel(r"$B_x$, nT", fontdict={"color": "k"})
        ax.set_ylim(-3000,3000)
        ax.axvline(frame.index.tolist()[-1000], ymin=4/6, ymax=5/6, color = "darkgreen", drawstyle="steps-mid")
        ax.axhline(2000, xmin=0.86, xmax=0.86+2e-2, color = "darkgreen")
        ax.axhline(1000, xmin=0.86, xmax=0.86+2e-2, color = "darkgreen")
        ax.text(frame.index.tolist()[-970], 1500, "1000 nT", ha="left", va="center", fontdict={"color": "darkgreen", "size":10})
        ax.set_yticklabels([])
        ax.legend(loc=2)
        ax = axes[1]
        ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
        ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
        ax.plot(frame.index, (base*multiplier[i])+frame.Y-np.mean(frame.Y), colors[i], ls="-", lw=1., label=stn.upper())
        ax.set_ylabel(r"$B_y$, nT", fontdict={"color": "k"})
        ax.axvline(frame.index.tolist()[2000], ymin=4/6, ymax=5/6, color = "darkgreen", drawstyle="steps-mid")
        ax.axhline(2000, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
        ax.axhline(1000, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
        ax.text(frame.index.tolist()[1970], 1500, "1000 nT", ha="right", va="center", fontdict={"color": "darkgreen", "size":10})
        ax.set_ylim(-3000,3000)
        ax.set_yticklabels([])
    axes[0].set_xlabel("Time, UT")
    axes[1].set_xlabel("Time, UT")
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig("docs/Bxy.Field.png", bbox_inches="tight")
    return

def plot_xy_electric_field_oneplot(stns, frames, dpi=150, wspace=0.2, hspace=0.1):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    fig, axes = plt.subplots(nrows=1, ncols=2, dpi=dpi, figsize=(12, 4))
    multiplier, colors = [4,3,2,1,0,-1,-2,-3,-4], ["r", "k", "b"]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
    base=100
    for i, frame in enumerate(frames):
        ax = axes[0]
        ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
        ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
        ax.plot(frame.index, (base*multiplier[i])+frame.X-np.mean(frame.X), colors[i], ls="-", lw=1.)
        ax.set_ylabel(r"$E_x$, mV/km", fontdict={"color": "k"})
        ax.set_ylim(-600,600)
        if i==0:
            ax.axvline(frame.index.tolist()[-1100], ymin=10/12, ymax=11/12, 
                       color = "darkgreen", drawstyle="steps-mid")
            ax.axhline(400, xmin=0.815, xmax=0.815+2e-2, color = "darkgreen")
            ax.axhline(500, xmin=0.815, xmax=0.815+2e-2, color = "darkgreen")
            ax.text(frame.index.tolist()[-1070], 450, "100 mv/km", ha="left", va="center", 
                    fontdict={"color": "darkgreen", "size":10})
        ax.set_yticklabels([])
        ax.legend(loc=2)
        ax = axes[1]
        ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
        ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
        txt = r"$Bin_{%d}[%s]$"%(i+1, stns[i].upper())
        ax.plot(frame.index, (base*multiplier[i])+frame.Y-np.mean(frame.Y), colors[i], ls="-", lw=1., label=txt)
        ax.set_ylabel(r"$E_y$, mV/km", fontdict={"color": "k"})
        if i==0:
            ax.axvline(frame.index.tolist()[2000], ymin=10/12, ymax=11/12, 
                       color = "darkgreen", drawstyle="steps-mid")
            ax.axhline(400, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
            ax.axhline(500, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
            ax.text(frame.index.tolist()[1970], 450, "100 mv/km", ha="right", va="center", 
                    fontdict={"color": "darkgreen", "size":10})
        #ax.axvline(frame.index.tolist()[2000], ymin=4/6, ymax=5/6, color = "darkgreen", drawstyle="steps-mid")
        #ax.axhline(2000, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
        #ax.axhline(1000, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
        #ax.text(frame.index.tolist()[1970], 1500, "1000 nT", ha="right", va="center", fontdict={"color": "darkgreen", "size":10})
        ax.set_ylim(-600,600)
        ax.set_yticklabels([])
    axes[0].set_xlabel("Time, UT")
    axes[1].set_xlabel("Time, UT")
    axes[1].legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig("docs/Exy.Field.png", bbox_inches="tight")
    return

def plot_xy_magnetic_field_fft(stns, frames, dpi=150, wspace=0.4, hspace=0.4):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    fig, axes = plt.subplots(nrows=len(stns), ncols=2, dpi=dpi, figsize=(12, 3*len(stns)))
    for i, stn in enumerate(stns):
        frame = frames[stn]
        ax = axes[i,0]
        ax.tick_params(axis="y", which="both", colors="b")
        ax.loglog(frame["X"]["freq"], np.absolute(frame["X"]["Xf"]), "b", ls="-", lw=1.)
        ax.text(0.05, 0.9, stn.upper()+"(x)", ha="left", va="center", transform=ax.transAxes, fontdict={"fontweight": "bold"})
        ax.set_ylabel(r"$|B_x(f)|$", fontdict={"color": "b"})
        ax.set_xlim(1e-6,1e-2)
        if i == len(stns)-1: ax.set_xlabel(r"Frequency, $f_0$")
        ax = ax.twinx()
        ax.tick_params(axis="y", which="both", colors="r")
        ax.spines["left"].set_color("b")
        ax.spines["right"].set_color("r")
        ax.semilogx(frame["X"]["freq"], np.rad2deg(np.angle(frame["X"]["Xf"])), "r", ls="-", lw=1., alpha=0.4)
        ax.set_ylim(-180, 180)
        ax.set_xlim(1e-6,1e-2)
        ax.set_ylabel(r"$\theta[B_x(f)]$", fontdict={"color": "r"})
        
        ax = axes[i,1]
        ax.tick_params(axis="y", which="both", colors="b")
        ax.loglog(frame["Y"]["freq"], np.absolute(frame["Y"]["Xf"]), "b", ls="-", lw=1.)
        ax.text(0.05, 0.9, stn.upper()+"(y)", ha="left", va="center", transform=ax.transAxes, fontdict={"fontweight": "bold"})
        ax.set_ylabel(r"$|B_y(f)|$", fontdict={"color": "b"})
        ax.set_xlim(1e-6,1e-2)
        if i == len(stns)-1: ax.set_xlabel(r"Frequency, $f_0$")
        ax = ax.twinx()
        ax.tick_params(axis="y", which="both", colors="r")
        ax.spines["left"].set_color("b")
        ax.spines["right"].set_color("r")
        ax.semilogx(frame["Y"]["freq"], np.rad2deg(np.angle(frame["Y"]["Xf"])), "r", ls="-", lw=1., alpha=0.4)
        ax.set_ylim(-180, 180)
        ax.set_xlim(1e-6,1e-2)
        ax.set_ylabel(r"$\theta[B_y(f)]$", fontdict={"color": "r"})
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig("docs/Bxy.Field.FFT.png", bbox_inches="tight")
    return

def plot_xy_electric_field_fft(stns, frames, dpi=150, wspace=0.4, hspace=0.4):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    fig, axes = plt.subplots(nrows=len(stns), ncols=2, dpi=dpi, figsize=(12, 3*len(stns)))
    for i, stn in enumerate(stns):
        frame = frames[stn]
        ax = axes[i,0]
        ax.tick_params(axis="y", which="both", colors="b")
        ax.loglog(frame["X"]["freq"], np.absolute(frame["X"]["Ef"]), "b", ls="-", lw=1.)
        ax.text(0.05, 0.9, stn.upper()+"(x)", ha="left", va="center", transform=ax.transAxes, fontdict={"fontweight": "bold"})
        ax.set_ylabel(r"$|E_x(f)|$", fontdict={"color": "b"})
        ax.set_xlim(1e-6,1e-2)
        if i == len(stns)-1: ax.set_xlabel(r"Frequency, $f_0$")
        ax = ax.twinx()
        ax.tick_params(axis="y", which="both", colors="r")
        ax.spines["left"].set_color("b")
        ax.spines["right"].set_color("r")
        ax.semilogx(frame["X"]["freq"], np.rad2deg(np.angle(frame["X"]["Ef"])), "r", ls="-", lw=1., alpha=0.4)
        ax.set_ylim(-180, 180)
        ax.set_xlim(1e-6,1e-2)
        ax.set_ylabel(r"$\theta[E_x(f)]$", fontdict={"color": "r"})
        
        ax = axes[i,1]
        ax.tick_params(axis="y", which="both", colors="b")
        ax.loglog(frame["Y"]["freq"], np.absolute(frame["Y"]["Ef"]), "b", ls="-", lw=1.)
        ax.text(0.05, 0.9, stn.upper()+"(y)", ha="left", va="center", transform=ax.transAxes, fontdict={"fontweight": "bold"})
        ax.set_ylabel(r"$|E_y(f)|$", fontdict={"color": "b"})
        ax.set_xlim(1e-6,1e-2)
        if i == len(stns)-1: ax.set_xlabel(r"Frequency, $f_0$")
        ax = ax.twinx()
        ax.tick_params(axis="y", which="both", colors="r")
        ax.spines["left"].set_color("b")
        ax.spines["right"].set_color("r")
        ax.semilogx(frame["Y"]["freq"], np.rad2deg(np.angle(frame["Y"]["Ef"])), "r", ls="-", lw=1., alpha=0.4)
        ax.set_ylim(-180, 180)
        ax.set_xlim(1e-6,1e-2)
        ax.set_ylabel(r"$\theta[E_y(f)]$", fontdict={"color": "r"})
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig("docs/Exy.Field.FFT.png", bbox_inches="tight")
    return


def plot_xy_electric_field(stns, frames, dpi=150, wspace=0.1, hspace=0.1):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    fig, axes = plt.subplots(nrows=2, ncols=1, dpi=dpi, figsize=(6, 6), 
                             sharex="all", sharey="all")
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
    for i, frame in enumerate(frames):
        col = colors[i]
        for j, lab in zip(range(2), ["x", "y"]):
            ax = axes[j]
            if i==0:
                ax.set_ylabel(r"$E_{%s}$, $mV.km^{-1}$"%lab)
                ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
                ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
                ax.xaxis.set_major_locator(mdates.DayLocator())
                ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
                ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 8)))
                #ax.set_ylim(-0.1,.2)
            txt = r"$Bin_{%d}[%s]$"%(i+1, stns[i].upper())
            ax.plot(frame.index, frame[lab.upper()], col, ls="-", lw=1., label=txt)
        if i == len(stns)-1: 
            ax.set_xlabel("Time, UT")
            ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig("docs/Exy.Field.png", bbox_inches="tight")
    return

def plot_induced_potential(stns, frames, dpi=150, wspace=0.1, hspace=0.1):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    fig, axes = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=(6, 3), 
                             sharex="all", sharey="all")
    ax = axes
    number = len(stns)
    cmap = plt.get_cmap("gnuplot")
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
    for i, frame in enumerate(frames):
        #ax = axes[i]
        col = "k"
        ax.set_ylabel(r"$V_j(t)=E_j^N(t)L_i^N+E_j^E(t)L_j^E$, $mV$", fontdict={"color": col})
        ax.spines["left"].set_color(col)
        ax.tick_params(axis="y", which="both", colors=col)
        ax.yaxis.label.set_color(col)
        ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
        ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 8)))
        
        #ax.set_ylim(-0.1,.2)
        txt = r"$Bin_{%d}[%s]$"%(i+1, stns[i].upper())
        ax.plot(frame.index, frame.Vj, color=colors[i], ls="-", lw=1., label=txt)
        if i == len(stns)-1: 
            ax.set_xlabel("Time, UT")
            ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig("docs/EField.png", bbox_inches="tight")
    return