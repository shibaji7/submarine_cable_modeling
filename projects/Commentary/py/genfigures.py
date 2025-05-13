import sys

sys.path.append("py/*")
import datetime as dt

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import utils

plt.style.use(["science", "ieee"])
plt.rcParams.update(
    {
        "text.usetex": False,
    }
)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]


def figure2():
    plt.rcParams.update({"font.size": 7})
    flim, M = [1e-10, 1e0], 1000
    freqs = np.linspace(flim[0], flim[1], int(flim[1] / (M * flim[0])) + 1)
    tfs = [
        utils.calcTFx(utils.sites[0], freqs),
        utils.calcTFx(utils.sites[1], freqs),
        utils.calcTF(utils.sites[2], freqs),
        utils.calcTF(utils.sites[3], freqs),
    ]

    fig = plt.figure(dpi=1000, figsize=(2.5, 3))
    ax = fig.add_subplot(211)
    ax.loglog(
        tfs[0].freq,
        np.abs(tfs[0].E2B),
        "r",
        lw=1.0,
        ls="-",
    )
    ax.loglog(tfs[2].freq, np.abs(tfs[2].E2B), "b", lw=0.6, ls="--", zorder=3)
    ax.loglog(tfs[3].freq, np.abs(tfs[3].E2B), "k", lw=0.6, ls="-", zorder=3)
    ax.set_ylabel(r"Amplitude [mV/km/nT]")
    ax.set_xlim(flim)
    ax.set_xticks([1e-10, 1e-5, 1e0])
    ax.set_ylim(1e-3, 1e3)
    ax.set_yticks([1e-3, 1e0, 1e3])
    ax.axvspan(1e-10, 1e-5, color="k", alpha=0.2)
    ax.set_xticks([])
    ax.text(
        0.9,
        0.9,
        rf"Case A, $\tau_1={utils.sites[0].layers[0].thickness/1000}$ km",
        ha="right",
        va="center",
        transform=ax.transAxes,
    )
    ax.text(0.1, 0.7, rf"(A-I)", ha="left", va="center", transform=ax.transAxes)
    ax = fig.add_subplot(212)
    ax.semilogx(
        tfs[0].freq, np.angle(tfs[0].E2B, deg=True), "r", lw=1.0, ls="-", label="non-U"
    )
    ax.semilogx(
        tfs[2].freq,
        np.angle(tfs[2].E2B, deg=True),
        "b",
        lw=0.6,
        ls="--",
        zorder=3,
        label=r"U($\rho$=3 $\Omega.m$)",
    )
    ax.semilogx(
        tfs[3].freq,
        np.angle(tfs[3].E2B, deg=True),
        "k",
        lw=0.6,
        ls="-",
        zorder=3,
        label=r"U($\rho$=3000 $\Omega.m$)",
    )
    ax.axvspan(1e-10, 1e-5, color="k", alpha=0.2)
    ax.set_yticks([0, 30, 45, 60, 90])
    ax.set_ylabel(r"Phase $[^\circ]$")
    ax.text(0.1, 0.7, rf"(A-II)", ha="left", va="center", transform=ax.transAxes)
    ax.legend(loc="upper right")
    ax.set_ylim(0, 90)
    ax.set_xticks([1e-10, 1e-5, 1e0])
    ax.set_xlim(flim)
    ax.set_xlabel(r"Frequency [Hz]")
    fig.subplots_adjust(wspace=0.15, hspace=0.15)
    fig.savefig("figures/Figure02.png", bbox_inches="tight")
    return


def figure3():
    plt.rcParams.update({"font.size": 7})
    flim, M = [1e-10, 1e0], 1000
    freqs = np.linspace(flim[0], flim[1], int(flim[1] / (M * flim[0])) + 1)
    tfs = [
        utils.calcTFx(utils.sites[0], freqs),
        utils.calcTFx(utils.sites[1], freqs),
        utils.calcTF(utils.sites[2], freqs),
        utils.calcTF(utils.sites[3], freqs),
    ]

    fig = plt.figure(dpi=1000, figsize=(2.5, 3))
    ax = fig.add_subplot(211)
    ax.loglog(tfs[1].freq, np.abs(tfs[1].E2B), "r", lw=1.0, ls="-")
    ax.loglog(tfs[2].freq, np.abs(tfs[2].E2B), "b", lw=0.6, ls="--", zorder=3)
    ax.loglog(tfs[3].freq, np.abs(tfs[3].E2B), "k", lw=0.6, ls="-", zorder=3)
    ax.axvspan(1e-10, 1e-5, color="k", alpha=0.2)
    ax.set_ylabel(r"Amplitude [mV/km/nT]")
    ax.text(0.1, 0.7, rf"(B-I)", ha="left", va="center", transform=ax.transAxes)
    ax.set_xlim(flim)
    ax.set_xticks([])
    ax.set_ylim(1e-3, 1e3)
    ax.set_yticks([1e-3, 1e0, 1e3])
    ax.text(
        0.9,
        0.9,
        rf"Case B, $\tau_1={utils.sites[1].layers[0].thickness/1000}$ km",
        ha="right",
        va="center",
        transform=ax.transAxes,
    )
    ax = fig.add_subplot(212)
    ax.semilogx(
        tfs[1].freq, np.angle(tfs[1].E2B, deg=True), "r", lw=1.0, ls="-", label="non-U"
    )
    ax.semilogx(
        tfs[2].freq,
        np.angle(tfs[2].E2B, deg=True),
        "b",
        lw=0.6,
        ls="--",
        zorder=3,
        label=r"U($\rho=3 \Omega.m$)",
    )
    ax.semilogx(
        tfs[3].freq,
        np.angle(tfs[3].E2B, deg=True),
        "k",
        lw=0.6,
        ls="-",
        zorder=3,
        label=r"U($\rho=3000 \Omega.m$)",
    )
    ax.axvspan(1e-10, 1e-5, color="k", alpha=0.2)
    ax.legend(loc="lower right")
    ax.set_yticks([0, 30, 45, 60, 90])
    ax.set_ylabel(r"Phase $[^\circ]$")
    ax.set_xlabel(r"Frequency [Hz]")
    ax.text(0.1, 0.7, rf"(B-II)", ha="left", va="center", transform=ax.transAxes)
    ax.set_ylim(0, 90)
    ax.set_xlim(flim)
    ax.set_xticks([1e-10, 1e-5, 1e0])
    fig.subplots_adjust(wspace=0.15, hspace=0.15)
    fig.savefig("figures/Figure03.png", bbox_inches="tight")
    return


def figure4():
    plt.rcParams.update({"font.size": 10})
    ds, _ = utils.get_benchmark_datasets()

    fig, start = plt.figure(dpi=1000, figsize=(5, 7)), 311
    ax = fig.add_subplot(start)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989, 3, 12, 12), dt.datetime(1989, 3, 15))
    ax.set_ylim(-4000, 2000)
    ax.set_ylabel(r"B, nT")
    ax.plot(ds.datetime, ds.x, ls="-", lw=0.6, color="r", label="$B_x$")
    ax.plot(ds.datetime, ds.y, ls="-", lw=0.6, color="k", label="$B_y$")
    ax.legend(loc=1)
    ax.set_xticklabels([])
    ax.text(0.05, 0.95, "(a)", ha="left", va="top", transform=ax.transAxes)

    # ax = fig.add_subplot(start+1)
    # date_format = mdates.DateFormatter("%H")
    # ax.xaxis.set_major_formatter(date_format)
    # ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    # ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    # ax.set_ylim(-500, 500)
    # # ax.set_ylabel(r"$\frac{\partial B}{\partial t}$ in frequency domain, nT/s")
    # ax.set_ylabel(r"GIC Index")
    # ax.plot(ds.datetime, utils.dBdt_fft(ds.x, sqrt=True), ls="-", lw=0.6, color="r")
    # ax.plot(ds.datetime, utils.dBdt_fft(ds.y, sqrt=True), ls="-", lw=0.6, color="k")
    # ax.text(0.05, 0.95, "(b)", ha="left", va="top", transform=ax.transAxes)
    # ax.set_xticklabels([])

    # ax = fig.add_subplot(start+2)
    # date_format = mdates.DateFormatter("%H")
    # ax.xaxis.set_major_formatter(date_format)
    # ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    # ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    # ax.set_ylim(-100, 100)
    # ax.set_ylabel(r"$\frac{\partial B}{\partial t}$ in frequency domain, nT/s")
    # # ax.set_ylabel(r"GIC Index")
    # ax.plot(datasets.datetime, dBdt_fft(datasets.x, sqrt=False)*a_scale, ls="-", lw=0.6, color="r")
    # ax.plot(datasets.datetime, dBdt_fft(datasets.y, sqrt=False)*a_scale, ls="-", lw=0.6, color="k")
    # ax.text(0.05, 0.95, "(c)", ha="left", va="top", transform=ax.transAxes)

    ax = fig.add_subplot(start + 1)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989, 3, 12, 12), dt.datetime(1989, 3, 15))
    ax.set_ylim(-30, 30)
    ax.set_ylabel(r"$\frac{\partial B}{\partial t}$ in time domain, nT/s")
    ax.set_xlabel(r"Time, UT since 12 UT on 12 March 1989")
    ax.plot(ds.datetime, ds.dx, ls="-", lw=0.6, color="r")
    ax.plot(ds.datetime, ds.dy, ls="-", lw=0.6, color="k")
    ax.text(0.05, 0.95, "(c)", ha="left", va="top", transform=ax.transAxes)
    fig.savefig("figures/Figure04.png", bbox_inches="tight")
    return


figure4()
figure2()
figure3()
