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


def plot_Bxy_stack(stns, frames, dpi=150, wspace=0.2, hspace=0.1):
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
    fig.savefig("prev/Bxy.Field.png", bbox_inches="tight")
    return

def plot_Exy_stack(stns, frames, dpi=150, wspace=0.2, hspace=0.1):
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
    fig.savefig("prev/Exy.png", bbox_inches="tight")
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