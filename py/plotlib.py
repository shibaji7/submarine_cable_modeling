"""plotlib.py: Module is used to implement plotting functions"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

import utility


def plot_Bxy_stack(stns, frames, dpi=150, wspace=0.2, hspace=0.1, fbase=""):
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
    fig.savefig(fbase + "Bxy.Field.png", bbox_inches="tight")
    return

def plot_Exy_stack(stns, frames, dpi=150, wspace=0.2, hspace=0.1, fbase=""):
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
        ax.set_ylim(-600,600)
        ax.set_yticklabels([])
    axes[0].set_xlabel("Time, UT")
    axes[1].set_xlabel("Time, UT")
    axes[1].legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig(fbase + "Exy.png", bbox_inches="tight")
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

def plot_induced_potential(stns, frames, dpi=150, wspace=0.1, hspace=0.1, fbase=""):
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
    fig.savefig(fbase + "EFieldx.png", bbox_inches="tight")
    return

def plot_total_potential(dx, dpi=150, wspace=0.1, hspace=0.1, fbase=""):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    fig, axes = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=(6, 3), 
                             sharex="all", sharey="all")
    ax = axes
    ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,12))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 8)))
    ax.plot(dx.index, dx["Vc(V)"], color="r", ls="-", lw=1., label=r"$\epsilon_c$")
    ax.plot(dx.index, dx["Vt(V)"], color="b", ls="-", lw=1., label=r"$V_c$")
    ax.set_xlabel("Time, UT")
    ax.set_ylabel("Voltage, V")
    ax.legend(loc=1)
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.savefig(fbase + "Pot.png", bbox_inches="tight")
    return

def potential_along_section(V, x, pname, sec=None, comp=None, Vi=None, Vk=None, Z=None, Y=None, gma=None, Z0=None, L=None):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    fig, axes = plt.subplots(nrows=1, ncols=1, dpi=150, figsize=(6, 3), 
                             sharex="all", sharey="all")
    ax = axes
    ax.set_ylabel("Voltage, V")
    ax.set_xlabel("Cable Length, km")
    ax.plot(x, V, "k", lw=0.8, ls="-")
    if (Z is not None): Z *= 1e3
    if (Y is not None): Y *= 1e3
    if (gma is not None): gma *= 1e3
    txt = ""
    if (sec is not None) and (comp is not None): txt += "Along: Bin%02d, %s\n"%(sec, str(comp))
    if (Vi is not None) and (Vk is not None): txt += (r"$V_i,V_k\sim %.1f V, %.1f V$"%(Vi,Vk) + "\n")
    if (Z is not None) and (Y is not None): txt += (r"$Z,Y\sim$ %s $\Omega/km$, %s $\mho/km$"%(utility.frexp102str(Z),utility.frexp102str(Y)) + "\n")
    if (gma is not None) and (Z0 is not None): txt += (r"$\gamma,Z_0\sim$ %s /km, %s $\Omega$"%(utility.frexp102str(gma),utility.frexp102str(Z0)) + "\n")
    if (L is not None): txt += "L=%d km"%L
    ax.text(0.05,0.95, txt, ha="left",va="top", transform=ax.transAxes, fontsize="small")
    ax.set_xlim(x[0], x[-1])
    fig.savefig(pname, bbox_inches="tight")
    return

def cable_potential(V, x, pname, comp="X"):
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    fig, axes = plt.subplots(nrows=1, ncols=1, dpi=150, figsize=(6, 3), 
                             sharex="all", sharey="all")
    ax = axes
    ax.set_ylabel("Voltage, V")
    ax.set_xlabel("Cable Length, km")
    ax.plot(x, V, "k", lw=0.8, ls="-")
    ax.set_xlim(x[0], x[-1])
    txt = ""
    if comp: txt += "Along: %s\n"%(str(comp))
    ax.text(0.05,0.95, txt, ha="left",va="top", transform=ax.transAxes, fontsize="small")
    fig.savefig(pname, bbox_inches="tight")
    return