"""
Diagnostic script: compare polar correlation plots under four conditions
to isolate the source of distortion in r(GIC, B_theta) and r(GIC, dB/dt).

2x2 grid per panel:
  top-left:     full time series, dtheta=1 deg
  top-right:    trimmed (10% each end), dtheta=1 deg
  bottom-left:  full time series, dtheta=10 deg
  bottom-right: trimmed (10% each end), dtheta=10 deg

Runs CaseA pipeline at 45 deg and CaseB pipeline at 45 deg.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scienceplots

sys.path.append("py/framework")

from intermagnet import Bfield
from mt_sites import PROFILES, TransferFunction
from pipeline import Pipeline
from network import Substation

plt.style.use(["science", "ieee"])
plt.rcParams.update({"text.usetex": False})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
mpl.rc("font", size=12)


def _corr_loop(bx, by, gic, dtheta, trim_fraction):
    """Core correlation loop (no class dependency)."""
    N = len(bx)
    i0 = int(N * trim_fraction)
    i1 = N - i0
    bx_t, by_t, gic_t = bx[i0:i1], by[i0:i1], gic[i0:i1]

    theta_deg = np.arange(0, 361, dtheta)
    cor = np.zeros(len(theta_deg))
    for i, th in enumerate(np.deg2rad(theta_deg)):
        bt = bx_t * np.cos(th) + by_t * np.sin(th)
        if np.std(bt) < 1e-10 or np.std(gic_t) < 1e-10:
            cor[i] = 0.0
        else:
            cor[i] = np.corrcoef(bt, gic_t)[0, 1]
    return theta_deg, cor


def _polar_ax(fig, pos, title="", color="r"):
    ax = fig.add_subplot(pos, polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rmax(1.0)
    ax.set_rticks([0, 0.5, 1.0])
    ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2])
    if title:
        ax.set_title(title, fontsize=9, pad=6)
    return ax


def make_2x2_panel(fig, subplot_spec, bx, by, gic, field_label, color, angle=45):
    """Fill a 2x2 block of polar axes in fig using the given subplot_spec."""
    from matplotlib.gridspec import GridSpecFromSubplotSpec

    inner = GridSpecFromSubplotSpec(2, 2, subplot_spec=subplot_spec,
                                    hspace=0.45, wspace=0.35)

    configs = [
        (0.0, 1.0, "full, 1°"),
        (0.1, 1.0, "trimmed, 1°"),
        (0.0, 10.0, "full, 10°"),
        (0.1, 10.0, "trimmed, 10°"),
    ]
    axes = []
    for idx, (trim, dth, label) in enumerate(configs):
        row, col = divmod(idx, 2)
        th, cor = _corr_loop(bx, by, gic, dth, trim)
        ax = fig.add_subplot(inner[row, col], polar=True)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.plot(np.deg2rad(th), np.abs(cor), color=color, lw=0.9)
        # Pipeline orientation line (solid) and perpendicular (dashed)
        if angle is not None:
            ax.plot(
                np.deg2rad([angle, angle + 180]), [1, 1],
                lw=1.2, ls="-", color="k",
            )
            ax.plot(
                np.deg2rad([angle + 90, angle + 270]), [1, 1],
                lw=0.6, ls="--", color="k",
            )
        ax.set_rmax(1.0)
        ax.set_rticks([0, 0.5, 1.0])
        ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2])
        ax.set_title(label, fontsize=8, pad=4)
        axes.append(ax)
    return axes


def run_pipeline_diagnostic(site_name, angle=45):
    site = getattr(PROFILES, site_name)
    bf = Bfield.create_benchmark_bfield()
    from mt_sites import TransferFunction
    tf = TransferFunction(site)
    ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta)
    pipe = Pipeline(angle=angle)
    gic = pipe.compute_J(pipe.compute_E(ef.ex, ef.ey))

    fig = plt.figure(figsize=(12, 10), dpi=150)
    fig.suptitle(
        f"{site_name} – Pipeline {angle}°: trimming & resolution comparison",
        fontsize=12, y=0.98,
    )

    from matplotlib.gridspec import GridSpec
    gs = GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.4)

    panel_configs = [
        (gs[0, 0], bf.bx, bf.by,  r"r(GIC, B$_\theta$)",        "r"),
        (gs[0, 1], bf.dbx, bf.dby, r"r(GIC, dB$_\theta$/dt)",    "b"),
    ]

    for spec, bx_arr, by_arr, label, clr in panel_configs:
        axes = make_2x2_panel(fig, spec, bx_arr, by_arr, gic, label, clr, angle=angle)
        # label the 2x2 block
        ax00 = axes[0]
        ax00.text(-0.25, 1.35, label, transform=ax00.transAxes,
                  fontsize=10, color=clr, fontweight="bold")

    outfile = f"figures/diag_polar_{site_name}_{angle}.png"
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outfile}")


def run_network_diagnostic(site_name, stn_names=None):
    """
    Generate the same 2x2 trimming/resolution diagnostic as the pipeline version,
    but for each substation in the power network.

    One output PNG per substation, e.g.:
        figures/diag_polar_net_CaseA_S2.png
    """
    if stn_names is None:
        stn_names = ["S2", "S3", "S4", "S5", "S6", "S8"]

    site = getattr(PROFILES, site_name)
    bf = Bfield.create_benchmark_bfield()
    tf = TransferFunction(site)
    ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta)

    from matplotlib.gridspec import GridSpec

    for stn in stn_names:
        net = Substation(name=stn)
        gic = net.compute_J(ef.ex, ef.ey)

        fig = plt.figure(figsize=(12, 10), dpi=150)
        fig.suptitle(
            f"{site_name} – Substation {stn}: trimming & resolution comparison",
            fontsize=12, y=0.98,
        )

        gs = GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.4)

        panel_configs = [
            (gs[0, 0], bf.bx,  bf.by,  r"r(GIC, B$_\theta$)",       "r"),
            (gs[0, 1], bf.dbx, bf.dby, r"r(GIC, dB$_\theta$/dt)",   "b"),
        ]

        for spec, bx_arr, by_arr, label, clr in panel_configs:
            # No pipeline angle line for substations — pass angle=None
            axes = make_2x2_panel(
                fig, spec, bx_arr, by_arr, gic, label, clr, angle=None
            )
            ax00 = axes[0]
            ax00.text(-0.25, 1.35, label, transform=ax00.transAxes,
                      fontsize=10, color=clr, fontweight="bold")

        outfile = f"figures/diag_polar_net_{site_name}_{stn}.png"
        fig.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {outfile}")


if __name__ == "__main__":
    for case in ["CaseA", "CaseB"]:
        run_pipeline_diagnostic(case, angle=45)
        run_network_diagnostic(case)
