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
    pot_time, params, omni,
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
    dpi=300, figsize=(6, 6), text_size=10,
    formatter = DateFormatter(r"%H^{%M}"),
):
    setups(text_size)
    fig, axes = plt.subplots(
        nrows=3, ncols=1, dpi=dpi, 
        figsize=figsize, sharex=True
    )
    ax = axes[0]
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_locator(major_locator)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.set_ylabel("B [nT]", fontdict={"color": "b"})
    ax.set_ylim(0, 100)
    ax.plot(
        omni.time, omni.B, 
        color="b", ls="-", 
        lw=linewdiths[0], alpha=0.7,
    )
    ax.set_xlim(dates)
    ax.text(
        0.05, 1.05, 
        fig_title, 
        ha="left", va="bottom", 
        transform=ax.transAxes
    )
    ax = ax.twinx()
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_locator(major_locator)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.set_ylabel(r"AE [nT]", fontdict={"color": "r"})
    ax.plot(
        omni.time, omni.AE, 
        color="r", ls="-", 
        lw=linewdiths[0], alpha=0.7,
    )
    ax.set_ylim(0, 2000)
    ax.set_xlim(dates)


    ax = axes[1]
    ax.xaxis.set_major_formatter(formatter)
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
    ax.legend(loc=3, prop={"size": 10})

    ax = axes[2]
    ax.xaxis.set_major_formatter(formatter)
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
    ax.set_xlim(xlim)
    ax.set_ylim(ylims[1])
    fig.savefig(fname, bbox_inches="tight")
    return


def plot_correlation_function(
    dct, X_params, Y_params,
    xlim, ylabels, ylim=[0, 1],
    fname="figures/correlation.png",
    xlabel="Time [UT]",
    dpi=300, text_size=10,
):
    setups(text_size)
    fig, axes = plt.subplots(
        nrows=len(X_params), ncols=1, dpi=dpi, 
        figsize=(8, len(X_params)*3),
    )
    major_locator=mdates.MinuteLocator(byminute=range(0, 60, 30))
    minor_locator=mdates.MinuteLocator(byminute=range(0, 60, 10))
    formatter = DateFormatter(r"%H^{%M}")
    for i, ax in enumerate(axes):
        if i < len(X_params)-1:
            ax.set_xticks([])
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_locator(major_locator)
        ax.xaxis.set_minor_locator(minor_locator)
        ax.set_ylabel(ylabels[i])
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    N, colors = len(Y_params), ["r", "b", "g", "k"]
    ax = axes[0]
    for i in range(N):
        ax.plot(dct["time"], dct[f"{X_params[0]}_{Y_params[i]}"], colors[i]+".", ms=0.6)

    ax = axes[1]
    for i in range(N):
        ax.plot(dct["time"], dct[f"{X_params[1]}_{Y_params[i]}"], colors[i]+".", ms=0.6)
    
    ax = axes[2]
    for i in range(N):
        ax.plot(dct["time"], dct[f"{X_params[2]}_{Y_params[i]}"], colors[i]+".", ms=0.6)
    
    ax = axes[3]
    for i in range(N):
        ax.plot(dct["time"], dct[f"{X_params[3]}_{Y_params[i]}"], colors[i]+".", ms=0.6)
    ax.set_xlabel(xlabel)
    
    fig.savefig(fname, bbox_inches="tight")
    return

def create_correlation_function(
    df, X_params, ylabels,
    fname="figures/correlation.png",
    dpi=300, text_size=12,
):
    setups(text_size)
    fig, axes = plt.subplots(
        nrows=2, ncols=2, dpi=dpi, 
        figsize=(8, 8), sharex="all"
    )
    axes = axes.ravel()

    ax = axes[0]
    ax.text(
        0.05, 1.05, "Shock @17:20 UT, on 10 May, 2024", 
        ha="left", va="center", transform=ax.transAxes
    )
    ax.plot(df["V(v)"], df[X_params[0]], "r.", ms=0.6, alpha=0.6)
    ax.set_ylabel(ylabels[0])
    ax = axes[1]
    ax.plot(df["V(v)"], df[X_params[1]], "r.", ms=0.6, alpha=0.6)
    ax.set_ylabel(ylabels[1])
    ax = axes[2]
    ax.plot(df["V(v)"], df[X_params[2]], "r.", ms=0.6, alpha=0.6)
    ax.set_ylabel(ylabels[2])
    ax = axes[3]
    ax.plot(df["V(v)"], df[X_params[3]], "r.", ms=0.6, alpha=0.6)
    ax.set_ylabel(ylabels[3])
    
    for i, ax in enumerate(axes):
        ax.set_xlim([-250, 150])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("Cable Voltage, V")

    fig.savefig(fname, bbox_inches="tight")
    return