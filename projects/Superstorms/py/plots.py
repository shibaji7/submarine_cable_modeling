import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt



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

class TimeSeriesPlot(object):

    def __init__(
        self, dates, fig_title, num_subplots=3, text_size=15,
        major_locator=mdates.HourLocator(byhour=range(0, 24, 12)),
        minor_locator=mdates.HourLocator(byhour=range(0, 24, 1)),
        formatter=DateFormatter(r"$%H^{%M}$"),
    ):
        self.dates = dates
        self.num_subplots = num_subplots
        self._num_subplots_created = 0
        self.text_size = text_size
        self.fig_title = fig_title
        self.major_locator = major_locator
        self.minor_locator = minor_locator
        self.formatter = formatter
        self.fig = plt.figure(figsize=(8, 3*num_subplots), dpi=180) # Size for website
        return

    def add_vlines(
        self, ax, vlines=[], colors=[]
    ):
        for col, vline in zip(colors, vlines):
            ax.axvline(vline, ls="--", lw=1.5, color=col)
            ax.axvline(vline, ls="--", lw=1.5, color=col)
        return
    
    def _add_axis(self):
        setups(self.text_size)
        self._num_subplots_created += 1
        ax = self.fig.add_subplot(self.num_subplots, 1, self._num_subplots_created)
        if self._num_subplots_created == 1:
            ax.text(
                0.075, 1.05, self.fig_title, 
                ha="left", va="center", 
                fontdict=dict(size=self.text_size),
                transform=ax.transAxes
            )
        ax.xaxis.set_major_locator(self.major_locator)
        ax.xaxis.set_minor_locator(self.minor_locator)
        ax.xaxis.set_major_formatter(self.formatter)
        return ax

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight")

    def close(self):
        self.fig.clf()
        plt.close()

    def add_omni(
        self, omni, ylelft="B", yright="FlowPressure", 
        linewdiths=[0.7, 1.5],
        ytitles = [r"B [nT]", r"$P_{dyn}$ [nPa]"],
        colors=["b", "r"], ylims=[[], []], xlabel=""
    ):
        ax = self._add_axis()
        ax.plot(
            omni.time, omni[ylelft], 
            color=colors[0], ls="-", 
            lw=linewdiths[0], alpha=0.7,
        )
        ax.set_ylabel(ytitles[0], fontdict=dict(color=colors[0]))
        ax.set_ylim(ylims[0])
        ax.set_xlabel(xlabel)
        ax.set_xlim(self.dates)

        axt =  ax.twinx()
        axt.xaxis.set_major_formatter(self.formatter)
        axt.xaxis.set_major_locator(self.major_locator)
        axt.xaxis.set_minor_locator(self.minor_locator)
        axt.plot(
            omni.time, omni[yright], 
            color=colors[1], ls="-", 
            lw=linewdiths[1], alpha=0.7,
        )
        axt.set_ylabel(ytitles[1], fontdict=dict(color=colors[1]))
        axt.set_ylim(ylims[1])
        return ax, axt

    def add_themis(
        self, themis_fgm, themis_mom, pnames=[],
        lw=0.7, colors=["k", "b", "g"],
        ylim=[-100, 100], xlabel="", loc=2,
        labels=[["$B_x$","$B_y$","$B_z$"]],
        ylabels=[r"$B_{sw}$ [nT]", r"$P_{dyn}$ [nPa]"]
    ):
        ax = self._add_axis()
        for j in range(themis_fgm[pnames[0]]["y"].shape[1]):
            ax.plot(
                [dt.datetime.utcfromtimestamp(x) for x in themis_fgm[pnames[0]]["x"]],
                themis_fgm[pnames[0]]["y"][:, j], color=colors[j], lw=lw, ls="-", label=labels[0][j]
            )
        ax.set_xlim(self.dates)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabels[0])
        ax.legend(loc=loc, prop={"size": 10})

        axt = ax.twinx()
        axt.xaxis.set_major_formatter(self.formatter)
        axt.xaxis.set_major_locator(self.major_locator)
        axt.xaxis.set_minor_locator(self.minor_locator)
        axt.plot(
            [dt.datetime.utcfromtimestamp(x) for x in themis_mom[pnames[1]]["x"]],
            themis_mom[pnames[1]]["y"], 
            color="m", ls="-", 
            lw=lw, alpha=0.7,
        )
        axt.set_xlim(self.dates)
        axt.set_ylabel(ylabels[1], fontdict=dict(color="m"))
        return ax
    
    def add_mag(
        self, df, stations=["FRD", "STJ", "HAD"],
        lw=0.7, colors=["k", "b", "g"],
        ylim=[-2000, 1500], xlabel="", loc=1,
        ylabel=r"$B_{x,y}$ [nT]"
    ):
        ax = self._add_axis()
        for c, stn in zip(colors, stations):
            o = df[stn]
            o.X = o.X - np.nanmean(o.X.iloc[:60*10])
            o.Y = o.Y - np.nanmean(o.Y.iloc[:60*10])
            o.Z = o.Z - np.nanmean(o.Z.iloc[:60*10])
            ax.plot(
                o.index, o.X, 
                color=c, ls="-", 
                lw=lw, alpha=0.7,
                label=fr"$B[{stn}]$"
            )
            ax.plot(
                o.index, o.Y, 
                color=c, ls="--", 
                lw=lw, alpha=0.7,
            )
        ax.set_xlim(self.dates)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(loc=loc, prop={"size": 10})
        return ax
    
    def add_voltage(
        self, df, lw=0.7, color="k",
        ylim=[-500, 500], xlabel="Time [UT]"
    ):
        ax = self._add_axis()
        ax.plot(
            df.index[::3], df["V(v)"][::3], 
            color=color, ls="-", 
            lw=lw, alpha=0.7,
        )
        ax.set_ylabel("Voltage, [V]")
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_xlim(self.dates)
        return ax

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
