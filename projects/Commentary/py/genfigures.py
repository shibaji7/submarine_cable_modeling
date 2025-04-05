import sys
sys.path.append("py/*")
import utils

import numpy as np

import matplotlib as mpl
import matplotlib.dates as mdates
import datetime as dt
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(["science", "ieee"])
plt.rcParams.update({
    "text.usetex": False,
})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]

def figure2():
    plt.rcParams.update({
        "font.size": 7
    })
    flim, M = [1e-10, 1e0], 1000
    freqs = np.linspace(flim[0], flim[1], int(flim[1] / (M*flim[0])) + 1)
    tfs = [
        utils.calcTFx(utils.sites[0], freqs),
        utils.calcTFx(utils.sites[1], freqs),
        utils.calcTF(utils.sites[2], freqs),
        utils.calcTF(utils.sites[3], freqs),
    ]


    fig = plt.figure(dpi=300, figsize=(2.5, 3))
    ax = fig.add_subplot(211)
    ax.loglog(tfs[0].freq, np.abs(tfs[0].E2B), "r", lw=1.0, ls="-",)
    ax.loglog(tfs[2].freq, np.abs(tfs[2].E2B), "b", lw=0.6, ls="--", zorder=3)
    ax.loglog(tfs[3].freq, np.abs(tfs[3].E2B), "k", lw=0.6, ls="-", zorder=3)
    ax.set_ylabel(r"Amplitude [mV/km/nT]")
    ax.set_xlim(flim)
    ax.set_xticks([1e-10, 1e-5, 1e0])
    ax.set_ylim(1e-3,1e3)
    ax.set_yticks([1e-3, 1e0, 1e3])
    ax.axvspan(1e-10, 1e-5, color="k", alpha=0.2)
    ax.set_xticks([])
    ax.text(0.9, 0.9, fr"Case A, $\tau_1={utils.sites[0].layers[0].thickness/1000}$ km", ha="right", va="center", transform=ax.transAxes)
    ax.text(0.1, 0.7, fr"(A-I)", ha="left", va="center", transform=ax.transAxes)
    ax = fig.add_subplot(212)
    ax.semilogx(tfs[0].freq, np.angle(tfs[0].E2B, deg=True), "r", lw=1.0, ls="-", label="non-U")
    ax.semilogx(tfs[2].freq, np.angle(tfs[2].E2B, deg=True), "b", lw=0.6, ls="--", zorder=3, label=r"U($\rho$=3 $\Omega.m$)")
    ax.semilogx(tfs[3].freq, np.angle(tfs[3].E2B, deg=True), "k", lw=0.6, ls="-", zorder=3, label=r"U($\rho$=3000 $\Omega.m$)")
    ax.axvspan(1e-10, 1e-5, color="k", alpha=0.2)
    ax.set_yticks([0, 30, 45, 60, 90])
    ax.set_ylabel(r"Phase $[^\circ]$")
    ax.text(0.1, 0.7, fr"(A-II)", ha="left", va="center", transform=ax.transAxes)
    ax.legend(loc='upper right')
    ax.set_ylim(0, 90)
    ax.set_xticks([1e-10, 1e-5, 1e0])
    ax.set_xlim(flim)
    ax.set_xlabel(r"Frequency [Hz]")
    fig.subplots_adjust(wspace=0.15,hspace=0.15)
    fig.savefig("figures/Figure02.png", bbox_inches="tight")
    return

def figure3():
    plt.rcParams.update({
        "font.size": 7
    })
    flim, M = [1e-10, 1e0], 1000
    freqs = np.linspace(flim[0], flim[1], int(flim[1] / (M*flim[0])) + 1)
    tfs = [
        utils.calcTFx(utils.sites[0], freqs),
        utils.calcTFx(utils.sites[1], freqs),
        utils.calcTF(utils.sites[2], freqs),
        utils.calcTF(utils.sites[3], freqs),
    ]

    fig = plt.figure(dpi=300, figsize=(2.5, 3))
    ax = fig.add_subplot(211)
    ax.loglog(tfs[1].freq, np.abs(tfs[1].E2B), "r", lw=1.0, ls="-")
    ax.loglog(tfs[2].freq, np.abs(tfs[2].E2B), "b", lw=0.6, ls="--", zorder=3)
    ax.loglog(tfs[3].freq, np.abs(tfs[3].E2B), "k", lw=0.6, ls="-", zorder=3)
    ax.axvspan(1e-10, 1e-5, color="k", alpha=0.2)
    ax.set_ylabel(r"Amplitude [mV/km/nT]")
    ax.text(0.1, 0.7, fr"(B-I)", ha="left", va="center", transform=ax.transAxes)
    ax.set_xlim(flim)
    ax.set_xticks([])
    ax.set_ylim(1e-3,1e3)
    ax.set_yticks([1e-3, 1e0, 1e3])
    ax.text(0.9, 0.9, fr"Case B, $\tau_1={utils.sites[1].layers[0].thickness/1000}$ km", ha="right", va="center", transform=ax.transAxes)
    ax = fig.add_subplot(212)
    ax.semilogx(tfs[1].freq, np.angle(tfs[1].E2B, deg=True), "r", lw=1., ls="-", label="non-U")
    ax.semilogx(tfs[2].freq, np.angle(tfs[2].E2B, deg=True), "b", lw=0.6, ls="--", zorder=3, label=r"U($\rho=3 \Omega.m$)")
    ax.semilogx(tfs[3].freq, np.angle(tfs[3].E2B, deg=True), "k", lw=0.6, ls="-", zorder=3, label=r"U($\rho=3000 \Omega.m$)")
    ax.axvspan(1e-10, 1e-5, color="k", alpha=0.2)
    ax.legend(loc='lower right')
    ax.set_yticks([0, 30, 45, 60, 90])
    ax.set_ylabel(r"Phase $[^\circ]$")
    ax.set_xlabel(r"Frequency [Hz]")
    ax.text(0.1, 0.7, fr"(B-II)", ha="left", va="center", transform=ax.transAxes)
    ax.set_ylim(0, 90)
    ax.set_xlim(flim)
    ax.set_xticks([1e-10, 1e-5, 1e0])
    fig.subplots_adjust(wspace=0.15,hspace=0.15)
    fig.savefig("figures/Figure03.png", bbox_inches="tight")
    return

def figure4():
    plt.rcParams.update({
        "font.size": 10
    })
    ds = utils.get_benchmark_datasets()

    fig, start = plt.figure(dpi=300, figsize=(5,7)), 311
    ax = fig.add_subplot(start)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    ax.set_ylim(-4000, 2000)
    ax.set_ylabel(r"B, nT")
    ax.plot(ds.datetime, ds.x, ls="-", lw=0.6, color="r", label="$B_x$")
    ax.plot(ds.datetime, ds.y, ls="-", lw=0.6, color="k", label="$B_y$")
    ax.legend(loc=1)
    ax.set_xticklabels([])
    ax.text(0.05, 0.95, "(a)", ha="left", va="top", transform=ax.transAxes)

    ax = fig.add_subplot(start+1)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    ax.set_ylim(-500, 500)
    # ax.set_ylabel(r"$\frac{\partial B}{\partial t}$ in frequency domain, nT/s")
    ax.set_ylabel(r"GIC Index")
    ax.plot(ds.datetime, utils.dBdt_fft(ds.x, sqrt=True), ls="-", lw=0.6, color="r")
    ax.plot(ds.datetime, utils.dBdt_fft(ds.y, sqrt=True), ls="-", lw=0.6, color="k")
    ax.text(0.05, 0.95, "(b)", ha="left", va="top", transform=ax.transAxes)
    ax.set_xticklabels([])

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

    ax = fig.add_subplot(start+2)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    ax.set_ylim(-300, 300)
    ax.set_ylabel(r"$\frac{\partial B}{\partial t}$ in time domain, nT/s")
    ax.set_xlabel(r"Time, UT since 12 UT on 12 March 1989")
    ax.plot(ds.datetime, ds.dx, ls="-", lw=0.6, color="r")
    ax.plot(ds.datetime, ds.dy, ls="-", lw=0.6, color="k")
    ax.text(0.05, 0.95, "(c)", ha="left", va="top", transform=ax.transAxes)
    fig.savefig("figures/Figure04.png", bbox_inches="tight")
    return

def find_BM_Efield():
    plt.rcParams.update({
        "font.size": 10
    })
    ds = utils.get_benchmark_datasets()
    site = utils.sites[4]

    from scubas.utils import fft, ifft
    Bxf, f = fft(ds.x, 1)
    Ets = dict(
        Y=ifft(np.array(utils.calcTFx(site, freqs=f).E2B) * Bxf),
    )
    Byf, f = fft(ds.y, 1)
    Ets["X"] = ifft(np.array(utils.calcTFx(site, freqs=f).E2B) * Byf)

    fig, start = plt.figure(dpi=300, figsize=(5,7)), 311
    ax = fig.add_subplot(start)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    ax.set_ylim(-4000, 2000)
    ax.set_ylabel(r"B, nT")
    ax.plot(ds.datetime, ds.x, ls="-", lw=0.6, color="r", label="$B_x$")
    ax.plot(ds.datetime, ds.y, ls="-", lw=0.6, color="k", label="$B_y$")
    ax.legend(loc=1)
    ax.set_xticklabels([])
    ax.text(0.05, 0.95, "(a)", ha="left", va="top", transform=ax.transAxes)

    ax = fig.add_subplot(start+1)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    ax.set_ylim(-10, 10)
    ax.set_ylabel(r"E, V/km""")
    ax.plot(ds.datetime, Ets["X"]/1e3, ls="-", lw=0.6, color="r")
    ax.plot(ds.datetime, Ets["Y"]/1e3, ls="-", lw=0.6, color="k")
    # ax.plot(ds.datetime, utils.dBdt_fft(ds.y, sqrt=True), ls="-", lw=0.6, color="k")
    ax.text(0.05, 0.95, "(b)", ha="left", va="top", transform=ax.transAxes)
    ax.set_xlabel(r"Time, UT since 12 UT on 12 March 1989")
    fig.savefig("figures/Figure04a.png", bbox_inches="tight")
    return

def figure5():
    from fan import CartoDataOverlay
    import cartopy
    cb = CartoDataOverlay(date=dt.datetime(2024,5,11), extent=[-90, -78, 28, 36],)
    ax = cb.add_axes()

    ax.add_ground_station("Sub1", 33.6135, -87.3737, markerColor="m", font_color="k", xOffset=-0.5, yOffset=0.)
    ax.add_ground_station("Sub2", 34.3104, -86.3658, markerColor="m", font_color="k", xOffset=-0.5, yOffset=0.25)
    ax.add_ground_station("Sub3", 33.9551, -84.6794, markerColor="m", font_color="k", xOffset=-0.5, yOffset=0.5)
    ax.add_ground_station("Sub4", 33.5479, -86.0746, markerColor="m", font_color="k", xOffset=-0.25, yOffset=-0.25)
    ax.add_ground_station("Sub5", 32.7051, -84.6634, markerColor="m", font_color="k", xOffset=-0.5, yOffset=-0.5)
    ax.add_ground_station("Sub6", 33.3773, -82.6188, markerColor="m", font_color="k", xOffset=-0.5, yOffset=-0.5)
    ax.add_ground_station("Sub7", 34.2522, -82.8363, markerColor="m", font_color="k", xOffset=-0.5, yOffset=0.5)
    ax.add_ground_station("Sub8", 34.1956, -81.0980, markerColor="m", font_color="k", xOffset=-0.5, yOffset=0.5)

    ax.add_bus(
        lats=[33.6135, 34.3104], lons=[-87.3737, -86.3658], color="b",
    )
    ax.add_bus(
        lats=[33.6135, 33.5479], lons=[-87.3737, -86.0746], color="b",
    )
    ax.add_bus(
        lats=[34.3104, 33.9551], lons=[-86.3658, -84.6794], color="b",
    )
    ax.add_bus(
        lats=[33.9551, 32.7051], lons=[-84.6794, -84.6634], color="b",
    )
    ax.add_bus(
        lats=[34.3104, 32.7051], lons=[-86.3658, -84.6634], color="b",
    )


    ax.add_bus(
        lats=[34.2522, 34.1956], lons=[-82.8363, -81.0980], color="r",
    )
    ax.add_bus(
        lats=[34.2522, 33.3773], lons=[-82.8363, -82.6188], color="r",
    )
    ax.add_bus(
        lats=[33.5479, 33.3773], lons=[-86.0746, -82.6188], color="r",
    )
    ax.add_bus(
        lats=[33.9551, 33.5479], lons=[-84.6794, -86.0746], color="r",
    )
    ax.add_bus(
        lats=[32.7051, 33.3773], lons=[-84.6634, -82.6188], color="r",
    )
    ax.add_bus(
        lats=[32.7051, 34.2522], lons=[-84.6634, -82.8363], color="r",
    )
    ax.add_bus(
        lats=np.array([33.3773, 33.9551])-5e-2, lons=np.array([-82.6188, -84.6794])-5e-2, color="r", lw=0.4
    )
    ax.add_bus(
        lats=np.array([33.3773, 33.9551])+5e-2, lons=np.array([-82.6188, -84.6794])+5e-2, color="r", lw=0.4
    )
    ax.add_bus(
        lats=np.array([33.5479, 32.7051])-5e-2, lons=np.array([-86.0746, -84.6634])-5e-2, color="r", lw=0.4
    )
    ax.add_bus(
        lats=np.array([33.5479, 32.7051])+5e-2, lons=np.array([-86.0746, -84.6634])+5e-2, color="r", lw=0.4
    )
    cb.save("figures/Figure05a.png")
    cb.close()
    return

def figure6():
    plt.rcParams.update({
        "font.size": 10
    })
    ds = utils.get_benchmark_datasets()
    # Case A
    site = utils.sites[0]
    from scubas.utils import fft, ifft
    Bxf, f = fft(ds.x, 1)
    Ets = dict(
        Y=ifft(np.array(utils.calcTFx(site, freqs=f).E2B) * Bxf),
    )
    Byf, f = fft(ds.y, 1)
    Ets["X"] = ifft(np.array(utils.calcTFx(site, freqs=f).E2B) * Byf)

    fig, start = plt.figure(dpi=300, figsize=(5,7)), 311
    ax = fig.add_subplot(start)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    ax.set_ylim(-4000, 2000)
    ax.set_ylabel(r"B, nT")
    ax.plot(ds.datetime, ds.x, ls="-", lw=0.6, color="r", label="$B_x$")
    ax.plot(ds.datetime, ds.y, ls="-", lw=0.6, color="k", label="$B_y$")
    ax.legend(loc=1)
    ax.set_xticklabels([])
    ax.text(0.05, 0.95, "(a)", ha="left", va="top", transform=ax.transAxes)

    ax = fig.add_subplot(start+1)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    ax.set_ylim(-2, 2)
    ax.set_ylabel(r"E, V/km")
    ax.plot(ds.datetime, Ets["X"]/1e3, ls="-", lw=0.6, color="r")
    ax.plot(ds.datetime, Ets["Y"]/1e3, ls="-", lw=0.6, color="k")
    ax.text(0.05, 0.95, "(b) "+r"Case A, $\tau_0=10.0$ km", ha="left", va="top", transform=ax.transAxes)
    r0, r1 = np.corrcoef(ds.x, Ets["Y"])[0,1], np.corrcoef(ds.y, Ets["X"])[0,1]
    ax.text(0.05, 0.3, "$r(B,E)_{x}$=%.3f"%r0, ha="left", va="top", transform=ax.transAxes, fontdict=dict(color="r",size=7))
    ax.text(0.05, 0.2, "$r(B,E)_{y}$=%.3f"%r1, ha="left", va="top", transform=ax.transAxes, fontdict=dict(color="k",size=7))
    r0, r1 = np.corrcoef(ds.dx, Ets["Y"])[0,1], np.corrcoef(ds.dy, Ets["X"])[0,1]
    ax.text(0.95, 0.3, "$r(\partial B,E)_{x}$=%.3f"%r0, ha="right", va="top", transform=ax.transAxes, fontdict=dict(color="r",size=7))
    ax.text(0.95, 0.2, "$r(\partial B,E)_{y}$=%.3f"%r1, ha="right", va="top", transform=ax.transAxes, fontdict=dict(color="k",size=7))
    ax.set_xticklabels([])

    ax = fig.add_subplot(start+2)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    ax.set_ylabel(r"GIC [Amps/$^\circ$]", color="b")
    ax.plot(ds.datetime, Ets["X"]*20.3*1e-3, ls="-", lw=0.3, color="b") # N, A/p, to bus 9
    ax.plot(ds.datetime, Ets["X"]*-27.67*1e-3, ls="--", lw=0.3, color="k") # N, A/p, to transformer T8
    ax.set_ylim([-60, 60])
    ax.set_xlabel(r"Time, UT since 12 UT on 12 March 1989")
    # ax = ax.twinx()
    # ax.set_ylabel(r"GIC [Amps]", color="r")
    # ax.plot(ds.datetime, Ets["X"]*-279.08*1e-3, ls="-", lw=0.3, color="r", alpha=0.5) # N, A, to substation S5
    # # ax.set_xlim([0, 24])
    # ax.set_ylim([-400, 400])
    ax.text(0.05, 0.95, "(c) GIC expecetd through B9 and T8 due to $E_x$", ha="left", va="top", transform=ax.transAxes)

    
    fig.savefig("figures/Figure06.png", bbox_inches="tight")
    return

def figure7():
    plt.rcParams.update({
        "font.size": 10
    })
    ds = utils.get_benchmark_datasets()
    # Case A
    site = utils.sites[1]
    from scubas.utils import fft, ifft
    Bxf, f = fft(ds.x, 1)
    Ets = dict(
        Y=ifft(np.array(utils.calcTFx(site, freqs=f).E2B) * Bxf),
    )
    Byf, f = fft(ds.y, 1)
    Ets["X"] = ifft(np.array(utils.calcTFx(site, freqs=f).E2B) * Byf)

    fig, start = plt.figure(dpi=300, figsize=(5,7)), 311
    ax = fig.add_subplot(start)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    ax.set_ylim(-4000, 2000)
    ax.set_ylabel(r"B, nT")
    ax.plot(ds.datetime, ds.x, ls="-", lw=0.6, color="r", label="$B_x$")
    ax.plot(ds.datetime, ds.y, ls="-", lw=0.6, color="k", label="$B_y$")
    ax.legend(loc=1)
    ax.set_xticklabels([])
    ax.text(0.05, 0.95, "(a)", ha="left", va="top", transform=ax.transAxes)

    ax = fig.add_subplot(start+1)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    ax.set_ylim(-30, 30)
    ax.set_ylabel(r"E, V/km")
    ax.plot(ds.datetime, Ets["X"]/1e3, ls="-", lw=0.6, color="r")
    ax.plot(ds.datetime, Ets["Y"]/1e3, ls="-", lw=0.6, color="k")
    ax.text(0.05, 0.95, "(b) "+r"Case B, $\tau_0=100.0$ km", ha="left", va="top", transform=ax.transAxes)
    r0, r1 = np.corrcoef(ds.x, Ets["Y"])[0,1], np.corrcoef(ds.y, Ets["X"])[0,1]
    ax.text(0.05, 0.3, "$r(B,E)_{x}$=%.3f"%r0, ha="left", va="top", transform=ax.transAxes, fontdict=dict(color="r",size=7))
    ax.text(0.05, 0.2, "$r(B,E)_{y}$=%.3f"%r1, ha="left", va="top", transform=ax.transAxes, fontdict=dict(color="k",size=7))
    r0, r1 = np.corrcoef(ds.dx, Ets["Y"])[0,1], np.corrcoef(ds.dy, Ets["X"])[0,1]
    ax.text(0.95, 0.3, "$r(\partial B,E)_{x}$=%.3f"%r0, ha="right", va="top", transform=ax.transAxes, fontdict=dict(color="r",size=7))
    ax.text(0.95, 0.2, "$r(\partial B,E)_{y}$=%.3f"%r1, ha="right", va="top", transform=ax.transAxes, fontdict=dict(color="k",size=7))
    ax.set_xticklabels([])

    ax = fig.add_subplot(start+2)
    date_format = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.set_xlim(dt.datetime(1989,3,12,12), dt.datetime(1989,3,15))
    ax.set_ylabel(r"GIC [Amps/$^\circ$]")
    ax.plot(ds.datetime, Ets["X"]*20.3*1e-3, ls="-", lw=0.3, color="b") # N, A/p, to bus 9
    ax.plot(ds.datetime, Ets["X"]*-27.67*1e-3, ls="--", lw=0.3, color="r", alpha=0.6) # N, A/p, to transformer T8
    ax.set_ylim([-500, 500])
    # ax.text(0.1, 0.9, tag[0], ha="left", va="center", transform=ax.transAxes)
    # ax.set_xlabel(r"Time, UT since 12 UT on 12 March 1989")
    # ax = ax.twinx()
    # ax.set_ylabel(r"GIC [Amps]", color="r")
    # ax.plot(ds.datetime, Ets["X"]*-279.08*1e-3, ls="-", lw=0.3, color="r", alpha=0.5) # N, A, to substation S5
    # # ax.set_xlim([0, 24])
    # ax.set_ylim([-4000, 4000])
    ax.text(0.05, 0.95, "(c) GIC expecetd through B9 and T8 due to $E_x$", ha="left", va="top", transform=ax.transAxes)
    fig.savefig("figures/Figure07.png", bbox_inches="tight")
    return

# figure7()
# figure6()
figure5()
# find_BM_Efield()
# figure4()
# figure2()
# figure3()
