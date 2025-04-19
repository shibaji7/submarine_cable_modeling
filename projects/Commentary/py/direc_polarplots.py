import sys
sys.path.append("py/*")
import utils

import numpy as np
from scubas.utils import fft, ifft
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
plt.rcParams.update({
    "font.size": 10
})

def smooth_filling(x, y, s):
    from scipy.interpolate import UnivariateSpline
    spl = UnivariateSpline(x, y, s=s)
    x_new = np.arange(361)
    y_new = spl(x_new)
    y_new[-1] = y_new[0]
    return x_new, spl(x_new)

def compute(snum=0):
    dates = [dt.datetime(1989,3,12,12), dt.datetime(1989,3,15)]
    ds, _, del_t = utils.run_benchmark_sim(snum=snum, dates=dates)
    ds["Eh"] = np.sqrt(ds.Ex**2+ds.Ey**2)
    ds["theta_Eh"] = np.degrees(np.pi+np.arctan2(ds.Ey, ds.Ex))
    ds["Bh"] = np.sqrt(ds.x_o**2+ds.y_o**2)
    ds["theta_Bh"] = np.degrees(np.pi+np.arctan2(ds.y_o, ds.x_o))
    ds["Bh"] = ds["Bh"] - np.mean(ds["Bh"].iloc[:60])
    ds["dBh"] = np.diff(ds["Bh"], prepend=ds["Bh"].iloc[0])/del_t
    ds["theta_dBh"] = np.degrees(np.pi+np.arctan2(ds.dy, ds.dx))
    return ds

substn = dict(
    stn2= dict(
        a=115.635, b=-189.290,
    ),
    stn3= dict(
        a=139.848, b=-109.492,
    ),
    stn4= dict(
        a=19.983, b=-124.582,
    ),
    stn5= dict(
        a=-279.077, b=-65.458,
    ),
    stn6= dict(
        a=-57.291, b=354.525,
    ),
    stn8= dict(
        a=60.902, b=134.298,
    ),
)


def set_last_mid_point(ser):
    a, b = ser[0], ser[-1]
    ser[0], ser[-1] = (a+b)/2, (a+b)/2
    return ser

def compute_segmented_corr(ds, dth, s0, s1, sub1, sub2):
    Z_pipe = 0.004921 # Ohm / km
    Cor_B, Cor_dB = [], []
    theta_bins = np.arange(dth, 360+dth, dth)
    tcs = []

    for tc in theta_bins:
        tcs.append(tc)
        tl, th = tc-dth, tc+dth
        o = ds[
            (ds.theta_Eh>=tl) & (ds.theta_Eh<th)
        ]
        E_A_pipe = np.cos(tc)*np.array(o.Eh)
        GIC_pipe = E_A_pipe/1e3/Z_pipe
        Cor_B.append(np.abs(np.corrcoef(np.abs(o.Bh), GIC_pipe)[0,1]))
        Cor_dB.append(np.abs(np.corrcoef(o.dBh, GIC_pipe)[0,1]))

    # Cor_B, Cor_dB = (
    #     set_last_mid_point(Cor_B), 
    #     set_last_mid_point(Cor_dB)
    # )
    # theta, Cor_B = smooth_filling(tcs, Cor_B, s0)
    # _, Cor_dB = smooth_filling(tcs, Cor_dB, s1)
    # Cor_B, Cor_dB = (
    #     set_last_mid_point(Cor_B), 
    #     set_last_mid_point(Cor_dB)
    # )
    # return tcs, Cor_B-sub1, Cor_dB-sub2
    return tcs, Cor_B, Cor_dB

def pipeline():
    dsA = compute(0)
    thetaA, CorA_B, CorA_dB = compute_segmented_corr(
        dsA, dth = 10, s0=0.035, 
        s1=0.15, sub1=3e-2, sub2=0
    )
    dsB = compute(1)
    thetaB, CorB_B, CorB_dB = compute_segmented_corr(
        dsB, dth = 10, s0=1e-4, 
        s1=0.1, sub1=0.4, sub2=5e-2
    )

    fig, axs = plt.subplots(2, 1, subplot_kw={"projection": "polar"}, dpi=300, figsize=(3,6))
    # axs[0].plot(np.deg2rad(thetaA), CorA_B, ls="-", color="r", label=r"r(B, GIC)")
    # axs[0].plot(np.deg2rad(thetaA), CorA_dB, ls="-", color="b", label=r"r($\partial$B, GIC)")
    axs[0].bar(np.deg2rad(thetaA), CorA_B, bottom=0.0, color="r", width=np.deg2rad(thetaA[1]-thetaA[0]), alpha=0.5, label=r"r(B, GIC)")
    axs[0].bar(np.deg2rad(thetaA), CorA_dB, bottom=0.0, color="b", width=np.deg2rad(thetaA[1]-thetaA[0]), alpha=0.5, label=r"r($\partial$B, GIC)")
    axs[0].set_rticks([0, 0.5, 1.0])
    axs[0].set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    axs[0].set_rmax(1)
    axs[0].set_rmin(0)
    axs[0].legend(bbox_to_anchor=(0.7, 1.2))
    axs[0].text(0.05, 1.2, f"Cases for Pipeline", ha="left", va="center", transform=axs[0].transAxes)
    axs[0].text(1.1, 0.5, "(a) "+r"Case A, $\tau_0=10.0$ km", ha="left", va="center", transform=axs[0].transAxes, rotation=90)

    # axs[1].plot(np.deg2rad(thetaB), CorB_B, ls="-", color="r", label=r"r(B, GIC)")
    # axs[1].plot(np.deg2rad(thetaB), CorB_dB, ls="-", color="b", label=r"r($\partial$B, GIC)")
    axs[1].bar(np.deg2rad(thetaB), CorB_B, bottom=0.0, color="r", width=np.deg2rad(thetaB[1]-thetaB[0]), alpha=0.5, label=r"r(B, GIC)")
    axs[1].bar(np.deg2rad(thetaB), CorB_dB, bottom=0.0, color="b", width=np.deg2rad(thetaB[1]-thetaB[0]), alpha=0.5, label=r"r($\partial$B, GIC)")
    
    axs[1].set_rticks([0, 0.5, 1.0])
    axs[1].set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    axs[1].set_rmax(1)
    axs[1].set_rmin(0)
    axs[1].text(1.1, 0.5, "(b) "+r"Case B, $\tau_0=100.0$ km", ha="left", va="center", transform=axs[1].transAxes, rotation=90)
    fig.savefig("figures/Figure11.png", bbox_inches="tight")
    return

def compute_segmented_stn_corr(stn, ds, dth, s0, s1, sub1, sub2):
    stnab = substn[stn]
    Cor_B, Cor_dB = [], []
    theta_bins = np.arange(dth, 360+dth, dth)
    tcs = []
    # ds["GIC"] = 1e-3*(stnab["a"]*np.cos(ds.theta_Eh)*np.array(ds.Eh)+\
    #         (stnab["b"]*np.sin(ds.theta_Eh)*np.array(ds.Eh)))
    ds["GIC"] = 1e-3*(stnab["a"]*np.array(ds.Ex)+\
            (stnab["b"]*np.array(ds.Ey)))
    
    fig, axs = plt.subplots(2, 1, dpi=300, figsize=(8,6))
    ax = axs[0]
    ax.plot(range(len(ds)), ds.Bh, lw=1.2, zorder=2)
    ax.set_ylabel("B, nT")
    ax = ax.twinx()
    ax.plot(range(len(ds)), ds.dBh, color="r", lw=0.2, zorder=1)
    ax.set_xlim(0, len(ds))
    ax.set_ylabel(r"$\partial$ B, nT/s", color="r")

    ax = axs[1]
    ax.set_xlim(0, len(ds))
    ax.plot(range(len(ds)), ds.Ex, lw=0.2, zorder=2, ls="--", color="green")
    ax.plot(range(len(ds)), ds.Ey, lw=0.2, zorder=2, ls=":", color="b")
    ax.plot(range(len(ds)), np.sqrt(ds.Ex**2+ds.Ey**2), lw=0.6, zorder=2, ls="-", color="k")
    ax.set_ylabel("E, mv/km")
    ax = ax.twinx()
    ax.plot(range(len(ds)), ds.GIC, lw=0.2, zorder=2, ls="--", color="r")
    ax.set_ylabel("GIC, A")
    
    fig.savefig("figures/Figure120.png", bbox_inches="tight")

    for tc in theta_bins:
        tcs.append(tc)
        tl, th = tc-dth, tc+dth
        o = ds[
            (ds.theta_Eh>=tl) & (ds.theta_Eh<th)
        ]
        # GIC = (stnab["a"]*np.cos(tc)*np.array(o.Eh)) +\
        #     (stnab["b"]*np.sin(tc)*np.array(o.Eh))
        Cor_B.append(np.abs(np.corrcoef(o.Bh, o.GIC)[0,1]))
        Cor_dB.append(np.abs(np.corrcoef(o.dBh, o.GIC)[0,1]))

    Cor_B, Cor_dB = (
        set_last_mid_point(Cor_B), 
        set_last_mid_point(Cor_dB)
    )
    # theta, Cor_B = smooth_filling(tcs, Cor_B, s0)
    # _, Cor_dB = smooth_filling(tcs, Cor_dB, s1)
    # Cor_B, Cor_dB = (
    #     set_last_mid_point(Cor_B), 
    #     set_last_mid_point(Cor_dB)
    # )
    # return tcs, Cor_B-sub1, Cor_dB-sub2
    return tcs, Cor_B, Cor_dB
    
def powernet(stn=6):
    dsA = compute(0)
    dsB = compute(1)

    thetaA, CorA_B, CorA_dB = compute_segmented_stn_corr(
        f"stn{stn}", dsA, dth = 5, s0=0.02, 
        s1=0.08, sub1=3e-2, sub2=0
    )

    thetaB, CorB_B, CorB_dB = compute_segmented_stn_corr(
        f"stn{stn}", dsB, dth = 5, s0=0.02, 
        s1=0.1, sub1=0.1, sub2=5e-2
    )
    
    fig, axs = plt.subplots(2, 1, subplot_kw={"projection": "polar"}, dpi=300, figsize=(3,6))
    # axs[0].plot(np.deg2rad(thetaA), CorA_B, ls="-", color="r", label=r"r(B, GIC)")
    axs[0].bar(np.deg2rad(thetaA), CorA_B, bottom=0.0, color="r", width=np.deg2rad(thetaA[1]-thetaA[0]), alpha=0.5, label=r"r(B, GIC)")
    axs[0].bar(np.deg2rad(thetaA), CorA_dB, bottom=0.0, color="b", width=np.deg2rad(thetaA[1]-thetaA[0]), alpha=0.5, label=r"r($\partial$B, GIC)")
    # axs[0].plot(np.deg2rad(thetaA), CorA_dB, ls="-", color="b", label=r"r($\partial$B, GIC)")
    axs[0].set_rticks([0, 0.5, 1.0])
    axs[0].set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    axs[0].set_rmax(1.)
    axs[0].set_rmin(0)
    axs[0].legend(bbox_to_anchor=(0.7, 1.2))
    axs[0].text(0.05, 1.2, f"Cases for Sub/{stn}", ha="left", va="center", transform=axs[0].transAxes)
    axs[0].text(1.1, 0.5, "(a) "+r"Case A, $\tau_0=10.0$ km", ha="left", va="center", transform=axs[0].transAxes, rotation=90)

    # axs[1].plot(np.deg2rad(thetaB), CorB_B, ls="-", color="r", label=r"r(B, GIC)")
    # axs[1].plot(np.deg2rad(thetaB), CorB_dB, ls="-", color="b", label=r"r($\partial$B, GIC)")
    axs[1].bar(np.deg2rad(thetaB), CorB_B, bottom=0.0, color="r", width=np.deg2rad(thetaB[1]-thetaB[0]), alpha=0.5, label=r"r(B, GIC)")
    axs[1].bar(np.deg2rad(thetaB), CorB_dB, bottom=0.0, color="b", width=np.deg2rad(thetaB[1]-thetaB[0]), alpha=0.5, label=r"r($\partial$B, GIC)")
    axs[1].set_rticks([0, 0.5, 1.0])
    axs[1].set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    axs[1].set_rmax(1.)
    axs[1].set_rmin(0)
    axs[1].text(1.1, 0.5, "(b) "+r"Case B, $\tau_0=100.0$ km", ha="left", va="center", transform=axs[1].transAxes, rotation=90)
    fig.savefig("figures/Figure12.png", bbox_inches="tight")
    return

pipeline()
powernet()