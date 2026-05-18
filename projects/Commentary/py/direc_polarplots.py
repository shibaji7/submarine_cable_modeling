"""
direc_polarplots.py — Legacy standalone polar plot script.

NOTE: This file uses its own inline correlation logic (not the framework classes).
The correlation computation has been rewritten to use proper B-field projections
B_theta(t) = Bx(t)*cos(theta) + By(t)*sin(theta) instead of binning by E-field
direction and correlating with scalar Bh.

If you are generating final paper figures, use pipeline_sim.py and network_sim.py
instead — they use the corrected framework classes.
"""

import sys

sys.path.append("py/*")
import datetime as dt

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
plt.rcParams.update({"font.size": 10})


def compute(snum=0):
    dates = [dt.datetime(1989, 3, 12, 12), dt.datetime(1989, 3, 15)]
    ds, _, del_t = utils.run_benchmark_sim(snum=snum, dates=dates)
    ds["Eh"] = np.sqrt(ds.Ex**2 + ds.Ey**2)
    ds["theta_Eh"] = np.degrees(np.pi + np.arctan2(ds.Ey, ds.Ex))
    # --- NOTE: Bh magnitude kept for reference, but not used for correlation ---
    ds["Bh"] = np.sqrt(ds.x_o**2 + ds.y_o**2)
    ds["Bh"] = ds["Bh"] - np.mean(ds["Bh"].iloc[:60])
    ds["dBh"] = np.diff(ds["Bh"], prepend=ds["Bh"].iloc[0]) / del_t
    return ds, del_t


substn = dict(
    stn2=dict(a=115.635, b=-189.290),
    stn3=dict(a=139.848, b=-109.492),
    stn4=dict(a=19.983, b=-124.582),
    stn5=dict(a=-279.077, b=-65.458),
    stn6=dict(a=-57.291, b=354.525),
    stn8=dict(a=60.902, b=134.298),
)


def compute_segmented_corr(ds, del_t, pipe_angle=45.0, dtheta=1.0):
    """
    Compute directional correlation between pipeline GIC and B_theta / dB_theta/dt.
    
    --- CHANGED: Complete rewrite ---
    OLD: Binned data by E-field direction, correlated GIC subsets with scalar Bh.
    NEW: For each azimuthal direction theta, project Bx,By onto theta, 
         compute correlation with full GIC time series.
    
    Parameters:
    ds: DataFrame with columns x (Bx), y (By), dx (dBx/dt), dy (dBy/dt), Ex, Ey
    del_t: time step in seconds
    pipe_angle: pipeline orientation angle in degrees from north
    dtheta: angular resolution in degrees
    
    Returns:
    theta_deg, Cor_B, Cor_dB: arrays of angles and correlation coefficients
    """
    Z_pipe = 0.004921  # Ohm/km

    # Compute pipeline GIC (fixed orientation)
    E_pipe = np.array(ds.Ex) * np.cos(np.deg2rad(pipe_angle)) + \
             np.array(ds.Ey) * np.sin(np.deg2rad(pipe_angle))
    GIC_pipe = E_pipe / 1e3 / Z_pipe

    # Get Bx, By and their time derivatives
    Bx = np.array(ds.x)  # tapered, mean-removed Bx
    By = np.array(ds.y)  # tapered, mean-removed By
    dBx = np.array(ds.dx)  # dBx/dt
    dBy = np.array(ds.dy)  # dBy/dt

    theta_deg = np.arange(0, 360, dtheta)
    Cor_B = np.zeros(len(theta_deg))
    Cor_dB = np.zeros(len(theta_deg))

    for i, th in enumerate(np.deg2rad(theta_deg)):
        # Project B and dB/dt onto direction theta
        B_theta = Bx * np.cos(th) + By * np.sin(th)
        dB_theta = dBx * np.cos(th) + dBy * np.sin(th)

        Cor_B[i] = np.corrcoef(B_theta, GIC_pipe)[0, 1]
        Cor_dB[i] = np.corrcoef(dB_theta, GIC_pipe)[0, 1]

    return theta_deg, Cor_B, Cor_dB


def compute_segmented_stn_corr(stn, ds, del_t, dtheta=1.0):
    """
    Compute directional correlation between substation GIC and B_theta / dB_theta/dt.
    
    --- CHANGED: Complete rewrite ---
    OLD: Binned data by E-field direction, correlated GIC subsets with scalar Bh.
    NEW: Same projection approach as pipeline version.
    """
    stnab = substn[stn]

    # Compute substation GIC
    GIC = 1e-3 * (stnab["a"] * np.array(ds.Ex) + stnab["b"] * np.array(ds.Ey))

    # Get Bx, By and their time derivatives
    Bx = np.array(ds.x)
    By = np.array(ds.y)
    dBx = np.array(ds.dx)
    dBy = np.array(ds.dy)

    theta_deg = np.arange(0, 360, dtheta)
    Cor_B = np.zeros(len(theta_deg))
    Cor_dB = np.zeros(len(theta_deg))

    for i, th in enumerate(np.deg2rad(theta_deg)):
        B_theta = Bx * np.cos(th) + By * np.sin(th)
        dB_theta = dBx * np.cos(th) + dBy * np.sin(th)

        Cor_B[i] = np.corrcoef(B_theta, GIC)[0, 1]
        Cor_dB[i] = np.corrcoef(dB_theta, GIC)[0, 1]

    return theta_deg, Cor_B, Cor_dB


def pipeline():
    dsA, del_tA = compute(0)
    # --- CHANGED: new function signature ---
    thetaA, CorA_B, CorA_dB = compute_segmented_corr(dsA, del_tA, pipe_angle=45.0)

    dsB, del_tB = compute(1)
    thetaB, CorB_B, CorB_dB = compute_segmented_corr(dsB, del_tB, pipe_angle=45.0)

    fig, axs = plt.subplots(
        2, 1, subplot_kw={"projection": "polar"}, dpi=300, figsize=(3, 6)
    )
    # --- CHANGED: use np.abs() for polar plotting ---
    axs[0].bar(
        np.deg2rad(thetaA),
        np.abs(CorA_B),
        bottom=0.0,
        color="r",
        width=np.deg2rad(thetaA[1] - thetaA[0]),
        alpha=0.5,
        label=r"r(B, GIC)",
    )
    axs[0].bar(
        np.deg2rad(thetaA),
        np.abs(CorA_dB),
        bottom=0.0,
        color="b",
        width=np.deg2rad(thetaA[1] - thetaA[0]),
        alpha=0.5,
        label=r"r($\partial$B, GIC)",
    )
    axs[0].set_rticks([0, 0.5, 1.0])
    axs[0].set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2])
    axs[0].set_rmax(1)
    axs[0].set_rmin(0)
    axs[0].legend(bbox_to_anchor=(1.2, 1.2))
    axs[0].text(
        0.05, 1.2, f"Cases for Pipeline",
        ha="left", va="center", transform=axs[0].transAxes,
    )
    axs[0].text(
        1.1, 0.5, "(a) " + r"Case A, $\tau_0=10.0$ km",
        ha="left", va="center", transform=axs[0].transAxes, rotation=90,
    )

    axs[1].bar(
        np.deg2rad(thetaB),
        np.abs(CorB_B),
        bottom=0.0,
        color="r",
        width=np.deg2rad(thetaB[1] - thetaB[0]),
        alpha=0.5,
        label=r"r(B, GIC)",
    )
    axs[1].bar(
        np.deg2rad(thetaB),
        np.abs(CorB_dB),
        bottom=0.0,
        color="b",
        width=np.deg2rad(thetaB[1] - thetaB[0]),
        alpha=0.5,
        label=r"r($\partial$B, GIC)",
    )
    axs[1].set_rticks([0, 0.5, 1.0])
    axs[1].set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2])
    axs[1].set_rmax(1)
    axs[1].set_rmin(0)
    axs[1].text(
        1.1, 0.5, "(b) " + r"Case B, $\tau_0=100.0$ km",
        ha="left", va="center", transform=axs[1].transAxes, rotation=90,
    )
    fig.savefig("figures/Figure11.png", bbox_inches="tight")
    return


def powernet(stn=6):
    dsA, del_tA = compute(0)
    dsB, del_tB = compute(1)

    # --- CHANGED: new function signature ---
    thetaA, CorA_B, CorA_dB = compute_segmented_stn_corr(f"stn{stn}", dsA, del_tA)
    thetaB, CorB_B, CorB_dB = compute_segmented_stn_corr(f"stn{stn}", dsB, del_tB)

    fig, axs = plt.subplots(
        2, 1, subplot_kw={"projection": "polar"}, dpi=300, figsize=(3, 6)
    )
    # --- CHANGED: use np.abs() for polar plotting ---
    axs[0].bar(
        np.deg2rad(thetaA),
        np.abs(CorA_B),
        bottom=0.0, color="r",
        width=np.deg2rad(thetaA[1] - thetaA[0]),
        alpha=0.5, label=r"r(B, GIC)",
    )
    axs[0].bar(
        np.deg2rad(thetaA),
        np.abs(CorA_dB),
        bottom=0.0, color="b",
        width=np.deg2rad(thetaA[1] - thetaA[0]),
        alpha=0.5, label=r"r($\partial$B, GIC)",
    )
    axs[0].set_rticks([0, 0.5, 1.0])
    axs[0].set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2])
    axs[0].set_rmax(1.0)
    axs[0].set_rmin(0)
    axs[0].legend(bbox_to_anchor=(0.7, 1.2))
    axs[0].text(
        0.05, 1.2, f"Cases for Sub/{stn}",
        ha="left", va="center", transform=axs[0].transAxes,
    )
    axs[0].text(
        1.1, 0.5, "(a) " + r"Case A, $\tau_0=10.0$ km",
        ha="left", va="center", transform=axs[0].transAxes, rotation=90,
    )

    axs[1].bar(
        np.deg2rad(thetaB),
        np.abs(CorB_B),
        bottom=0.0, color="r",
        width=np.deg2rad(thetaB[1] - thetaB[0]),
        alpha=0.5, label=r"r(B, GIC)",
    )
    axs[1].bar(
        np.deg2rad(thetaB),
        np.abs(CorB_dB),
        bottom=0.0, color="b",
        width=np.deg2rad(thetaB[1] - thetaB[0]),
        alpha=0.5, label=r"r($\partial$B, GIC)",
    )
    axs[1].set_rticks([0, 0.5, 1.0])
    axs[1].set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2])
    axs[1].set_rmax(1.0)
    axs[1].set_rmin(0)
    axs[1].text(
        1.1, 0.5, "(b) " + r"Case B, $\tau_0=100.0$ km",
        ha="left", va="center", transform=axs[1].transAxes, rotation=90,
    )
    fig.savefig("figures/Figure12.png", bbox_inches="tight")
    return


pipeline()
# powernet()
