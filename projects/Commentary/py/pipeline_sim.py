import datetime as dt
import sys

sys.path.append("py/framework")

from intermagnet import Bfield
from mt_sites import PROFILES, TransferFunction
from pipeline import Pipeline, sqrtf_weighted_signal
from plots import StackPlots

import numpy as np


def run_pipe_line_simulations_for_benchmark_event(
    site_name: str,
    title: str = "",
    angle: float = 90,
    date_lim: list[dt.datetime] = [
        dt.datetime(1989, 3, 12, 12),
        dt.datetime(1989, 3, 15),
    ],
):
    """
    Run the pipeline simulation for the benchmark event.
    """
    site = getattr(PROFILES, site_name)
    bf = Bfield.create_benchmark_bfield()
    tf = TransferFunction(site)
    ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta)
    pipeline = Pipeline(angle=angle)
    e_pipe = pipeline.compute_E(ef.ex, ef.ey)
    gic_pipe = pipeline.compute_J(e_pipe)

    sp = StackPlots(nrows=2, ncols=2, dpi=300)
    _, ax = sp.plot_stack_plots(
        bf.date,
        gic_pipe,
        dict(start=(0, 0), colspan=2, rowspan=1),
        title=rf"{title} / Pipe oriented: {angle}$^\circ$",
        ylabel="GIC (A)",
        color="k",
        lw=0.9,
        xlim=date_lim,
        ylim=[-1500, 1500],
        datetime=True,
        tag="(a)",
        xlabel="Time, UT"
    )

    # --- CHANGED: pass bx, by instead of bf.bmag ---
    theta, cor = pipeline.compute_segmented_correlation(
        bf.bx, bf.by, ef.ex, ef.ey
    )
    sp.plot_dirctional_plots(
        theta,
        np.abs(cor),
        dict(start=(1, 0), colspan=1, rowspan=1),
        text=r"r(GIC, $B_\theta$)",
        color="r",
        cable_angle=angle,
        tag="(b)",
    )

    # --- CHANGED: pass dbx, dby instead of bf.dbh ---
    theta, cor = pipeline.compute_segmented_correlation(
        bf.dbx, bf.dby, ef.ex, ef.ey
    )
    sp.plot_dirctional_plots(
        theta,
        np.abs(cor),
        dict(start=(1, 1), colspan=1, rowspan=1),
        text=r"r(GIC, $\frac{\partial B_\theta}{\partial t}$)",
        color="b",
        cable_angle=angle,
        tag="(c)",
    )
    sp.fig.subplots_adjust(hspace=0.5)
    sp.save_fig(f"figures/Comp_Paper_Benchmark_Pipeline_{site_name}_{angle}.png")
    sp.close()
    return


def _print_proxy_summary(angle, theta_deg, curves):
    """
    Print a numeric summary table for the proxy robustness test.

    Parameters
    ----------
    angle : float
        Pipeline orientation [degrees].
    theta_deg : np.ndarray
        Azimuth array.
    curves : list of (label, r_array) tuples
    """
    print(f"\nProxy robustness summary — Pipeline at {angle}°")
    print(f"{'Curve':<40} {'Peak |r|':>8}  {'Peak θ (°)':>10}")
    print("-" * 62)
    for label, r in curves:
        idx = np.argmax(r)
        print(f"{label:<40} {r[idx]:>8.4f}  {theta_deg[idx]:>10.1f}")
    print()


def run_proxy_robustness_pipeline(
    angle: float = 45,
    date_lim: list[dt.datetime] = [
        dt.datetime(1989, 3, 12, 12),
        dt.datetime(1989, 3, 15),
    ],
):
    """
    Reviewer 1, Comment 1.1(c) — sqrt(f) proxy robustness test for pipeline.

    Computes GIC under Case A and Case B from the same benchmark Bx/By, then
    correlates a single conductivity-agnostic sqrt(f)·B_theta proxy against
    both GIC series across all azimuths.  Generates one figure per case
    (matching filenames Comp_Paper_Benchmark_Pipeline_CaseA/B_<angle>.png)
    with the two proxy curves overlaid as dashed lines on the existing B_theta
    and dB/dt polar panels.  Also prints a numeric summary table.
    """
    bf = Bfield.create_benchmark_bfield()

    # --- Compute Case A and Case B E-fields and GIC ---
    ef_A = TransferFunction(getattr(PROFILES, "CaseA")).compute_Efield(
        bf.bx, bf.by, bf.del_ta
    )
    ef_B = TransferFunction(getattr(PROFILES, "CaseB")).compute_Efield(
        bf.bx, bf.by, bf.del_ta
    )
    pipe = Pipeline(angle=angle)
    gic_A = pipe.compute_J(pipe.compute_E(ef_A.ex, ef_A.ey))
    gic_B = pipe.compute_J(pipe.compute_E(ef_B.ex, ef_B.ey))

    # --- Proxy robustness (conductivity-agnostic) ---
    theta, r_proxy_A, r_proxy_B = pipe.compute_proxy_robustness(
        bf.bx, bf.by, gic_A, gic_B, bf.del_ta
    )

    # Proxy colours: m (vs GIC_A) and darkgreen (vs GIC_B)
    PROXY_A_COLOR = "m"   # m
    PROXY_B_COLOR = "darkgreen"   # darkgreen

    j = 0
    # --- Figure ---
    sp = StackPlots(nrows=2, ncols=2, dpi=300)
    for site_name, title, ef, gic_self, gic_other, self_label, other_label in [
        ("CaseA", "Case A", ef_A, gic_A, gic_B,
         r"r(GIC, $\sqrt{f}B_\theta$)", r"r(GIC, $\sqrt{f}B_\theta$)"),
        ("CaseB", "Case B", ef_B, gic_B, gic_A,
         r"r(GIC, $\sqrt{f}B_\theta$)", r"r(GIC, $\sqrt{f}B_\theta$)"),
    ]:
        gic_pipe = gic_A if site_name == "CaseA" else gic_B

        # --- Existing B_theta and dB/dt correlations for this case ---
        theta_B, cor_B = pipe.compute_segmented_correlation(
            bf.bx, bf.by, ef.ex, ef.ey
        )
        theta_dB, cor_dB = pipe.compute_segmented_correlation(
            bf.dbx, bf.dby, ef.ex, ef.ey
        )

        # --- Numeric summary ---
        # proxy curves are the same for both cases (same bx/by);
        # r_proxy_A = |r(proxy, GIC_A)|, r_proxy_B = |r(proxy, GIC_B)|
        r_self  = r_proxy_A if site_name == "CaseA" else r_proxy_B
        r_cross = r_proxy_B if site_name == "CaseA" else r_proxy_A
        _print_proxy_summary(
            angle,
            theta,
            [
                (f"r(GIC, B_theta)  [{site_name}]",  np.abs(cor_B)),
                (f"r(GIC, dB/dt)    [{site_name}]",  np.abs(cor_dB)),
                (self_label,                          r_self),
                (other_label,                         r_cross),
            ],
        )

        # B_theta panel with proxy overlaid
        _, ax_B = sp.plot_dirctional_plots(
            theta_B,
            np.abs(cor_B),
            dict(start=(j, 0), colspan=1, rowspan=1),
            text=r"r(GIC, $B_\theta$)" if j==0 else "",
            color="r",
            cable_angle=angle,
            tag="(a)",
        )
        sp.plot_dirctional_plots(
            theta_B, r_self,
            ax=ax_B,
            text=self_label if j==0 else "",
            grid=None,
            color=PROXY_A_COLOR,
            text_location=(1., 0.8),
            lw=0.6,
            tag="(b)",
        )

        # dB/dt panel with proxy overlaid
        _, ax_dB = sp.plot_dirctional_plots(
            theta_dB,
            np.abs(cor_dB),
            dict(start=(j, 1), colspan=1, rowspan=1),
            text=r"r(GIC, $\frac{\partial B_\theta}{\partial t}$)" if j==0 else "",
            color="b",
            cable_angle=angle,
            tag="(c)",
        )
        sp.plot_dirctional_plots(
            theta_B, r_cross,
            ax=ax_dB,
            grid=None,
            color=PROXY_A_COLOR,
            text_location=(0.6, 1.1),
            lw=0.6,
            tag="(d)",
        )
        if j==0:
            ax_B.text(0.5, 1.15, "Case A", ha="center", va="bottom", transform=ax_B.transAxes)
            ax_dB.text(0.5, 1.15, "Case B", ha="center", va="bottom", transform=ax_dB.transAxes)
        j += 1

    sp.fig.subplots_adjust(hspace=0.5, wspace=0.1)
    sp.save_fig(
        f"figures/Comp_Paper_Benchmark_Pipeline_{angle}.png"
    )
    sp.close()

    return


if __name__ == "__main__":
    for angle in [45]:
        run_pipe_line_simulations_for_benchmark_event(
            "CaseA", r"Case A", angle=angle
        )
        run_pipe_line_simulations_for_benchmark_event(
            "CaseB", r"Case B", angle=angle
        )
    run_proxy_robustness_pipeline(angle=angle)
