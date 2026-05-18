import datetime as dt
import sys

sys.path.append("py/framework")

from intermagnet import Bfield
from mt_sites import PROFILES, TransferFunction
from pipeline import Pipeline
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


if __name__ == "__main__":
    for angle in [45]:
        run_pipe_line_simulations_for_benchmark_event(
            "CaseA", r"Case A", angle=angle
        )
        run_pipe_line_simulations_for_benchmark_event(
            "CaseB", r"Case B", angle=angle
        )
