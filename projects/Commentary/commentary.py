import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import numpy as np
import pandas as pd
from scipy import constants as C

from scubas.datasets import Site
from scubas.models import OceanModel
from scubas.utils import fft, ifft


def calcZ(site, freqs):
    omega = 2 * C.pi * freqs
    sigma_s = 1 / site.layers[0].resistivity
    k2 = 1.0j * omega * C.mu_0 * sigma_s
    k = np.sqrt(k2)
    Zo = (1.0j * omega * C.mu_0 / k) / (C.mu_0 / 1.0e-3)
    return Zo

def calcTF(site, freqs):
    """
    Calculate the transfer functions.
    """
    Zo = calcZ(site, freqs)
    Zd = Zo
    omega = 2 * C.pi * freqs
    sigma_s = 1 / site.layers[0].resistivity
    k2 = 1.0j * omega * C.mu_0 * sigma_s
    k = np.sqrt(k2)
    kd = k * 0

    func = Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Zo))
    tf = pd.DataFrame()
    tf["freq"], tf["E2B"] = freqs, func
    return tf

def calcTFx(site, freqs):
    """
    Calculate the transfer functions.
    """
    Zo = calcZ(site, freqs)
    Zd = site.calcZ(freqs)[1]
    omega = 2 * C.pi * freqs
    sigma_s = 1 / site.layers[0].resistivity
    k2 = 1.0j * omega * C.mu_0 * sigma_s
    k = np.sqrt(k2)
    kd = k * 0

    func = Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Zo))
    tf = pd.DataFrame()
    tf["freq"], tf["E2B"] = freqs, func
    return tf

def abline():
    gca = plt.gca()
    gca.set_autoscale_on(False)
    gca.plot(gca.get_xlim(),gca.get_ylim(), ls="-.", lw=0.4, color="k")
    return

def abxline():
    gca = plt.gca()
    gca.set_autoscale_on(False)
    gca.plot(gca.get_xlim(),[4e-1,4e-1], ls="-.", lw=0.4, color="k")
    return

sites = [
    Site.init(
        conductivities=[1/3, 1/3000],
        thicknesses=[10000, np.inf],
        names=[
            "Sediments",
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper",
        site_name="Two layer Earth model",
    ),
    Site.init(
        conductivities=[1/3000, 1/3],
        thicknesses=[100000, np.inf],
        names=[
            "Sediments",
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper",
        site_name="Two layer Earth model",
    ),
    Site.init(
        conductivities=[1/3],
        thicknesses=[np.inf],
        names=[
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper",
        site_name="Unform Earth",
    ),
    Site.init(
        conductivities=[1/3000],
        thicknesses=[np.inf],
        names=[
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper",
        site_name="Unform Earth",
    ),
]


flim, M = [1e-10, 1e0], 1000
freqs = np.linspace(flim[0], flim[1], int(flim[1] / (M*flim[0])) + 1)
tfs = [
    calcTFx(sites[0], freqs),
    calcTFx(sites[1], freqs),
    calcTF(sites[2], freqs),
    calcTF(sites[3], freqs),
]

fig = plt.figure(dpi=300, figsize=(2.5, 3))
ax = fig.add_subplot(211)
ax.loglog(tfs[0].freq, np.abs(tfs[0].E2B), "r", lw=1.0, ls="-",)
ax.loglog(tfs[2].freq, np.abs(tfs[2].E2B), "b", lw=0.6, ls="--", zorder=3)
ax.loglog(tfs[3].freq, np.abs(tfs[3].E2B), "k", lw=0.6, ls="-", zorder=3)
#ax.legend(loc='upper left', bbox_to_anchor=(1.1, 1.0),)
#ax.set_xlabel(r"Frequency [Hz]")
ax.set_ylabel(r"Amplitude [mV/km/nT]")
ax.set_xlim(flim)
ax.set_xticks([1e-10, 1e-5, 1e0])
ax.set_ylim(1e-3,1e3)
ax.set_yticks([1e-3, 1e0, 1e3])
ax.axvspan(1e-10, 1e-5, color="k", alpha=0.2)
#ax.yaxis.tick_left()
ax.set_xticks([])
ax.text(0.9, 0.9, fr"Case A, $\tau_1={sites[0].layers[0].thickness/1000}$ km", ha="right", va="center", transform=ax.transAxes)
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
#ax.yaxis.set_label_position("right")
#ax.yaxis.tick_right()
ax.set_ylim(0, 90)
ax.set_xticks([1e-10, 1e-5, 1e0])
#ax.set_xticks([])
ax.set_xlim(flim)
ax.set_xlabel(r"Frequency [Hz]")
fig.subplots_adjust(wspace=0.15,hspace=0.15)
fig.savefig("cmnt1.png", bbox_inches="tight")

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
ax.text(0.9, 0.9, fr"Case B, $\tau_1={sites[1].layers[0].thickness/1000}$ km", ha="right", va="center", transform=ax.transAxes)
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
#ax.yaxis.set_label_position("right")
ax.set_ylim(0, 90)
ax.set_xlim(flim)
ax.set_xticks([1e-10, 1e-5, 1e0])
fig.subplots_adjust(wspace=0.15,hspace=0.15)
fig.savefig("cmnt11.png", bbox_inches="tight")


fig = plt.figure(dpi=300, figsize=(2.5, 3))
ax = fig.add_subplot(211)
ax.loglog(tfs[2].freq, np.abs(tfs[2].E2B), "b", lw=0.6, ls="--", zorder=3)
ax.loglog(tfs[3].freq, np.abs(tfs[3].E2B), "k", lw=1., ls="-", zorder=3)
ax.axvspan(1e-10, 1e-5, color="k", alpha=0.2)
ax.set_ylabel(r"Amplitude [mV/km/nT]")
ax.text(0.1, 0.7, fr"(a)", ha="left", va="center", transform=ax.transAxes)
ax.set_xlim(flim)
ax.set_xticks([])
ax.set_ylim(1e-4,1e2)
ax.set_yticks([1e-4, 1e-1, 1e2])
ax = fig.add_subplot(212)
ax.semilogx(tfs[2].freq, np.angle(tfs[2].E2B, deg=True), "b", lw=0.6, ls="--", zorder=3, label=r"U($\rho=3 \Omega.m$)")
ax.semilogx(tfs[3].freq, np.angle(tfs[3].E2B, deg=True), "k", lw=1., ls="-", zorder=1, label=r"U($\rho=3000 \Omega.m$)")
ax.set_yticks([0, 30, 45, 60, 90])
ax.axvspan(1e-10, 1e-5, color="k", alpha=0.2)
ax.legend(loc='lower right')
ax.set_ylabel(r"Phase $[^\circ]$")
ax.set_xlabel(r"Frequency [Hz]")
ax.text(0.1, 0.7, fr"(b)", ha="left", va="center", transform=ax.transAxes)
#ax.yaxis.set_label_position("right")
ax.set_ylim(0, 90)
ax.set_xlim(flim)
ax.set_xticks([1e-10, 1e-5, 1e0])
fig.subplots_adjust(wspace=0.15,hspace=0.15)
fig.savefig("cmnt12.png", bbox_inches="tight")


t = np.linspace(0, 24*3600, num=24*3600, endpoint=False)
B = np.zeros_like(t)
Am, Tm_min, Phim = (
    [ 
        200, 90, 30,
        17, 8, 3.5
    ], 
    np.array([ 
        180, 80, 36, 
        15, 8, 3
    ])*60,
    np.deg2rad([
        10, 20, 30, 40,
        50, 60
    ])
)
for i in range(len(Am)):
    B = B + ( Am[i] * np.sin((2*np.pi*t/Tm_min[i]) + Phim[i]) )
    
Bf, f = fft(B, 1)
Ets = [
    ifft(np.array(calcTFx(sites[0], freqs=f).E2B) * Bf),
    ifft(np.array(calcTFx(sites[1], freqs=f).E2B) * Bf),
    ifft(np.array(calcTFx(sites[2], freqs=f).E2B) * Bf),
    ifft(np.array(calcTFx(sites[3], freqs=f).E2B) * Bf)
]

tag = [
    r"Case A, $\tau_0=10.0$ km", r"Case B, $\tau_0=100.0$ km",
    r"U($\rho$=3 $\Omega.m$)", r"U($\rho$=3000 $\Omega.m$)"
]
fig = plt.figure(dpi=300, figsize=(4,1.5*3))
ax = fig.add_subplot(311)
ax.set_ylabel(r"B(t) [nT]", color="r")
ax.set_xlim([0, 24])
ax.set_ylim(-400, 400)
ax.set_xticks([0, 6, 12, 18, 24])
ax.set_xticklabels([])
ax.plot(t/3600, B, ls="-", lw=0.6, color="r", zorder=2)

ax = fig.add_subplot(312)
ax.set_ylabel(r"E(t) [mv/km]", color="b")
ax.plot(t/3600, Ets[0], ls="-", lw=0.6, color="b")
ax.set_xlim([0, 24])
ax.set_xticks([0, 6, 12, 18, 24])
r = np.corrcoef(B, Ets[0])[0,1]
ax.text(0.1, 0.9, tag[0], ha="left", va="center", transform=ax.transAxes)
ax.set_xticklabels([])
ax.text(0.9, 0.9, r"$r_{E,B}=%.3f$"%r, ha="right", va="center", transform=ax.transAxes)
ax.set_ylim(-200, 200)

ax = fig.add_subplot(313)
ax.set_ylabel(r"E(t) [mv/km]", color="b")
ax.plot(t/3600, Ets[1], ls="-", lw=0.6, color="b")
ax.set_xlim([0, 24])
ax.set_xticks([0, 6, 12, 18, 24])
r = np.corrcoef(B, Ets[1])[0,1]
ax.text(0.1, 0.9, tag[1], ha="left", va="center", transform=ax.transAxes)
ax.text(0.9, 0.9, r"$r_{E,B}=%.3f$"%r, ha="right", va="center", transform=ax.transAxes)
ax.set_ylim(-200, 200)
ax.set_xlabel("Time [UT]")
fig.subplots_adjust(wspace=0.05,hspace=0.05)
fig.savefig("cmnt2.png", bbox_inches="tight")



fig = plt.figure(dpi=300, figsize=(4,1.5*3))
dB = np.diff(B, prepend=B[0])
ax = fig.add_subplot(311)
ax.set_ylabel(r"$\frac{\partial B}{\partial t}$ [nT/s]", color="darkgreen")
ax.set_xlim([0, 24])
ax.set_ylim(-0.75, 0.75)
ax.set_xticks([0, 6, 12, 18, 24])
ax.set_xticklabels([])
ax.plot(t/3600, dB, ls="-", lw=0.6, color="darkgreen", zorder=2)

ax = fig.add_subplot(312)
ax.set_ylabel(r"E(t) [mv/km]", color="b")
ax.plot(t/3600, Ets[0], ls="-", lw=0.6, color="b")
ax.set_xlim([0, 24])
ax.set_xticks([0, 6, 12, 18, 24])
r = np.corrcoef(dB, Ets[0])[0,1]
ax.text(0.1, 0.9, tag[0], ha="left", va="center", transform=ax.transAxes)
ax.set_xticklabels([])
ax.text(0.9, 0.9, r"$r_{E,\partial B/\partial t}=%.3f$"%r, ha="right", va="center", transform=ax.transAxes)
ax.set_ylim(-200, 200)

ax = fig.add_subplot(313)
ax.set_ylabel(r"E(t) [mv/km]", color="b")
ax.plot(t/3600, Ets[1], ls="-", lw=0.6, color="b")
ax.set_xlim([0, 24])
ax.set_xticks([0, 6, 12, 18, 24])
r = np.corrcoef(dB, Ets[1])[0,1]
ax.text(0.1, 0.9, tag[1], ha="left", va="center", transform=ax.transAxes)
ax.text(0.9, 0.9, r"$r_{E,\partial B/\partial t}=%.3f$"%r, ha="right", va="center", transform=ax.transAxes)
ax.set_ylim(-200, 200)
ax.set_xlabel("Time [UT]")
fig.subplots_adjust(wspace=0.05,hspace=0.05)
fig.savefig("cmnt3.png", bbox_inches="tight")

fig = plt.figure(dpi=300, figsize=(4,1.5*2))
ax = fig.add_subplot(211)
ax.set_ylabel(r"GIC [Amps/$^\circ$]", color="b")
ax.plot(t/3600, Ets[0]*20.3*1e-3, ls="-", lw=0.6, color="b") # N, A/p, to bus 9
ax.plot(t/3600, Ets[0]*-27.67*1e-3, ls=":", lw=0.6, color="b") # N, A/p, to transformer T8
ax.set_xlim([0, 24])
ax.set_ylim([-6, 6])
ax.text(0.1, 0.9, tag[0], ha="left", va="center", transform=ax.transAxes)
ax.set_xticklabels([])
ax = ax.twinx()
ax.set_ylabel(r"GIC [Amps]", color="r")
ax.plot(t/3600, Ets[0]*-279.08*1e-3, ls="-", lw=0.3, color="r") # N, A, to substation S5
ax.set_xlim([0, 24])
ax.set_ylim([-40, 40])
ax.set_xticks([0, 6, 12, 18, 24])

ax = fig.add_subplot(212)
ax.set_ylabel(r"GIC [Amps/$^\circ$]", color="b")
ax.plot(t/3600, Ets[1]*20.3*1e-3, ls="-", lw=0.6, color="b") # N, A/p, to bus 9
ax.plot(t/3600, Ets[1]*-27.67*1e-3, ls=":", lw=0.6, color="b") # N, A/p, to transformer T8
ax.set_xlim([0, 24])
ax.set_ylim([-6, 6])
ax.text(0.1, 0.9, tag[1], ha="left", va="center", transform=ax.transAxes)
ax.set_xlabel("Time [UT]")
ax = ax.twinx()
ax.set_ylabel(r"GIC [Amps]", color="r")
ax.plot(t/3600, Ets[1]*-279.08*1e-3, ls="-", lw=0.3, color="r") # N, A, to substation S5
ax.set_xlim([0, 24])
# ax.set_ylim([-40, 40])
ax.set_xticks([0, 6, 12, 18, 24])
fig.subplots_adjust(wspace=0.05,hspace=0.05)
fig.savefig("GIC.png", bbox_inches="tight")


# fig = plt.figure(dpi=300, figsize=(4,1.5*3))
# dB = np.diff(B, prepend=B[0])
# ax = fig.add_subplot(311)
# ax.plot(t/3600, B, ls="-", lw=0.9, color="r", zorder=2)
# ax.set_ylabel(r"B(t) [nT]", color="r")
# ax.set_ylim(-400, 400)
# ax = ax.twinx()
# ax.set_ylabel(r"$\frac{\partial}{\partial t}$ B(t) [nT/s]", color="darkgreen")
# ax.set_xlim([0, 24])
# ax.set_ylim(-0.75, 0.75)
# ax.set_xticks([0, 6, 12, 18, 24])
# ax.set_xticklabels([])
# ax.plot(t/3600, dB, ls="-", lw=0.6, color="darkgreen", zorder=0)

# ax = fig.add_subplot(312)
# ax.set_ylabel(r"E(t) [mv/km]", color="b")
# ax.plot(t/3600, Ets[0], ls="-", lw=0.6, color="b")
# ax.set_xlim([0, 24])
# ax.set_xticks([0, 6, 12, 18, 24])
# ax.text(0.1, 0.9, tag[0], ha="left", va="center", transform=ax.transAxes)
# ax.set_xticklabels([])
# r = np.corrcoef(B, Ets[0])[0,1]
# ax.text(0.9, 0.9, r"$r_{E,B}=%.2f$"%r, ha="right", va="center", transform=ax.transAxes, fontdict={"color":"r"})
# r = np.corrcoef(dB, Ets[0])[0,1]
# ax.text(0.9, 0.8, r"$r_{E,\frac{\partial B}{\partial t}}=%.2f$"%r, ha="right", va="center", transform=ax.transAxes, fontdict={"color":"darkgreen"})
# ax.set_ylim(-300, 300)

# ax = fig.add_subplot(313)
# ax.set_ylabel(r"E(t) [mv/km]", color="b")
# ax.plot(t/3600, Ets[1], ls="-", lw=0.6, color="b")
# ax.set_xlim([0, 24])
# ax.set_xticks([0, 6, 12, 18, 24])
# ax.text(0.1, 0.9, tag[1], ha="left", va="center", transform=ax.transAxes)
# r = np.corrcoef(B, Ets[1])[0,1]
# ax.text(0.9, 0.9, r"$r_{E,B}=%.2f$"%r, ha="right", va="center", transform=ax.transAxes, fontdict={"color":"r"})
# r = np.corrcoef(dB, Ets[1])[0,1]
# ax.text(0.9, 0.8, r"$r_{E,\frac{\partial B}{\partial t}}=%.2f$"%r, ha="right", va="center", transform=ax.transAxes, fontdict={"color":"darkgreen"})
# ax.set_ylim(-300, 300)
# ax.set_xlabel("Time [UT]")
# fig.subplots_adjust(wspace=0.05,hspace=0.05)
# fig.savefig("cmnt3.png", bbox_inches="tight")