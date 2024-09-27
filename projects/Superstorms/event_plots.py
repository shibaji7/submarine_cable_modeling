import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def setups(size = 12):
    plt.style.use(["science", "ieee"])
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [
        "Tahoma", "DejaVu Sans",
        "Lucida Grande", "Verdana"
    ]
    mpl.rcParams.update(
        {"xtick.labelsize": size, "ytick.labelsize": size, "font.size": size}
    )
    return

def plot_potential(
    frames, dates, 
    pot_time, params,
    fname="figures/Pot.png",
    fig_title="Date: 10-12 May 2024", 
    linewdiths=[0.7, 1.5], alphas=[0.7],
    yticks=[[-2000, -1000, 0, 1000]],
    ylims=[[-2000, 1000], [-500, 500]],
    major_locator=mdates.HourLocator(byhour=range(0, 24, 12)),
    minor_locator=mdates.HourLocator(byhour=range(0, 24, 1)),
    ytitles = [r"$B_{x,y}$ [nT]", "Potential [Volts]"], 
    stations=["FRD", "STJ", "HAD"],
    colors=[["k", "b", "g"], ["k"]], xlabel="Time [UT]",
    dpi=300, figsize=(6, 4), text_size=10,
):
    setups(text_size)
    fig, axes = plt.subplots(
        nrows=2, ncols=1, dpi=dpi, 
        figsize=figsize, sharex=True
    )
    ax = axes[0]
    ax.xaxis.set_major_formatter(DateFormatter(r"%H^{%M}"))
    ax.xaxis.set_major_locator(major_locator)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.set_ylabel(ytitles[0])
    ax.set_ylim(ylims[0])
    for c, stn in zip(colors[0], stations):
        o = frames[stn]
        o.X = o.X - np.nanmean(o.X.iloc[:60*10])
        o.Y = o.Y - np.nanmean(o.Y.iloc[:60*10])
        o.Z = o.Z - np.nanmean(o.Z.iloc[:60*10])
        ax.plot(
            o.index, o.X, 
            color=c, ls="-", 
            lw=linewdiths[0], alpha=0.7,
            label=fr"$B[{stn}]$"
        )
        ax.plot(
            o.index, o.Y, 
            color=c, ls="--", 
            lw=linewdiths[0], alpha=0.7,
        )
    ax.set_yticks(yticks[0])
    ax.set_xlim(dates)
    ax.legend(loc=4)
    ax.text(
        0.05, 1.05, 
        fig_title, 
        ha="left", va="bottom", 
        transform=ax.transAxes
    )

    ax = axes[1]
    ax.xaxis.set_major_formatter(DateFormatter(r"%H^{%M}"))
    ax.xaxis.set_major_locator(major_locator)
    ax.xaxis.set_minor_locator(minor_locator)
    for c, i, alpha in zip(colors[1], range(len(params)), alphas):
        ax.plot(
            pot_time, 
            params[i], 
            color=c,
            ls="-", lw=linewdiths[1], 
            alpha=alpha
        )
    ax.set_ylabel(ytitles[1])
    ax.set_ylim(ylims[1])
    ax.set_xlim(dates)
    ax.set_xlabel(xlabel)
    fig.savefig(fname, bbox_inches="tight")
    return

def plot_transfer_functions(
    df, xlim=[1e-6,1e0], ylims=[[1e-3,1e3],[-90,90]],
    fname="figures/Transfer.png",
    ylabels = [r"Amplitude [mV/km/nT]", r"Phase $[^\circ]$"],
    yticks = [[1e-3, 1e-0, 1e3],[-90, -45, 0, 45, 90]],
    xlabel="Frequency [Hz]",
    xticks=[1e-6, 1e-3, 1e0],
    dpi=300, figsize=(4, 4), text_size=10,
):
    setups(text_size)
    fig, ax = plt.subplots(
        nrows=1, ncols=1, dpi=dpi, 
        figsize=figsize,
    )
    ax.loglog(df["stats"].freqs, df["stats"].amp, "k", lw=1.0, ls="-", alpha=1)
    ax.loglog(df["stats"].freqs, df["stats"].amp_ub_1, "k", lw=0.7, ls="--", alpha=0.7)
    ax.loglog(df["stats"].freqs, df["stats"].amp_lb_1, "k", lw=0.7, ls="--", alpha=0.7)
    ax.loglog(df["stats"].freqs, df["stats"].amp_ub_2, "k", lw=0.5, ls=":", alpha=0.4)
    ax.loglog(df["stats"].freqs, df["stats"].amp_lb_2, "k", lw=0.5, ls=":", alpha=0.4)
    ax.set_ylabel(ylabels[0])
    ax.set_xlim(xlim)
    ax.set_ylim(ylims[0])
    ax.set_xlabel(xlabel)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks[0])
    ax = ax.twinx()
    ax.semilogx(df["stats"].freqs, df["stats"].phase, "r", lw=1.0, ls="-", alpha=1)
    ax.semilogx(df["stats"].freqs, df["stats"].phase_ub_1, "r", lw=0.7, ls="--", alpha=0.7)
    ax.semilogx(df["stats"].freqs, df["stats"].phase_lb_1, "r", lw=0.7, ls="--", alpha=0.7)
    ax.semilogx(df["stats"].freqs, df["stats"].phase_ub_2, "r", lw=0.5, ls=":", alpha=0.4)
    ax.semilogx(df["stats"].freqs, df["stats"].phase_lb_2, "r", lw=0.5, ls=":", alpha=0.4)
    ax.set_ylabel(ylabels[1], fontdict=dict(color="r"))
    ax.set_ylim(ylims[1])
    fig.savefig(fname, bbox_inches="tight")
    return
