import datetime as dt
import sys

sys.path.append("py/framework")

from intermagnet import Bfield
from mt_sites import PROFILES, TransferFunction
from pipeline import Pipeline
from plots import StackPlots


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

    This function simulates the geoelectric field (E-field) and geomagnetically induced currents (GIC)
    along a pipeline for a benchmark geomagnetic event. It uses a predefined magnetic field (B-field)
    dataset, computes the corresponding E-field using a transfer function for the specified site,
    and calculates the resulting GIC along the pipeline. The results are visualized in a series of
    plots, including time series of the B-field, E-field, and GIC, as well as directional sensitivity
    plots showing correlations between GIC and magnetic field variations.

    Args:
        site_name (str): The site identifier/name for which the transfer function is applied.
        title (str): A descriptive title for the plots, indicating the site or scenario being analyzed.
        angle (float): The angle of the pipeline relative to the geomagnetic north (in degrees).
        date_lim (list[dt.datetime]): The time range for the simulation and plots.
    """
    site = getattr(PROFILES, site_name)
    bf = Bfield.create_benchmark_bfield()
    tf = TransferFunction(site)
    ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta)
    pipeline = Pipeline(angle=angle)
    e_pipe = pipeline.compute_E(ef.ex, ef.ey)
    gic_pipe = pipeline.compute_J(e_pipe)

    sp = StackPlots(nrows=2, ncols=2)
    # _, ax = sp.plot_stack_plots(
    #     bf.date,
    #     ef.ex / 1e3,
    #     dict(start=(0, 0), colspan=2, rowspan=1),
    #     text=f"E-field under {title}",
    #     ylabel="E-field (V/km)",
    #     color="r",
    #     lw=0.4,
    #     label=r"$E_x$",
    #     tag="(A)",
    #     xlim=date_lim,
    #     ylim=[-2, 2],
    #     datetime=True,
    # )
    # sp.plot_stack_plots(
    #     bf.date,
    #     ef.ey / 1e3,
    #     None,
    #     label=r"$E_y$",
    #     color="k",
    #     lw=0.4,
    #     xlim=date_lim,
    #     ax=ax,
    #     ylim=[-2, 2],
    #     datetime=True,
    # )
    # ax.legend(loc=1)
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
        tag="(A)",
        xlabel="Time, UT"
    )
    # sp.plot_stack_plots(
    #     bf.date,
    #     gic_pipe,
    #     None,
    #     ylabel="GIC (A)",
    #     color="r",
    #     lw=0.4,
    #     xlim=date_lim,
    #     ax=ax.twinx(),
    #     ylabel_color="r",
    #     ylim=[-500, 500],
    #     datetime=True,
    # )

    theta, cor = pipeline.compute_segmented_correlation(
        bf.bmag, ef.ex, ef.ey, normalize=site_name == "CaseA"
    )
    sp.plot_dirctional_plots(
        theta,
        cor,
        dict(start=(1, 0), colspan=1, rowspan=1),
        # title=rf"{title} / Pipe oriented: {angle}$^\circ$",
        text=r"r(GIC, $B_h$)",
        color="r",
        cable_angle=angle,
        tag="(B)",
    )
    theta, cor = pipeline.compute_segmented_correlation(
        bf.dbh, ef.ex, ef.ey, normalize=site_name == "CaseB"
    )
    sp.plot_dirctional_plots(
        theta,
        cor,
        dict(start=(1, 1), colspan=1, rowspan=1),
        text=r"r(GIC, $\partial B_h$)",
        color="b",
        cable_angle=angle,
        tag="(C)",
    )
    sp.fig.subplots_adjust(hspace=0.5)
    sp.save_fig(f"figures/Comp_Paper_Benchmark_Pipeline_{site_name}_{angle}.png")
    sp.close()
    return


if __name__ == "__main__":
    # for angle in np.arange(0, 100, 15):
    for angle in [45]:
        run_pipe_line_simulations_for_benchmark_event(
            "CaseA", r"Case A, $\tau_1$=10.0 km", angle=angle
        )
        run_pipe_line_simulations_for_benchmark_event(
            "CaseB", r"Case B, $\tau_1$=100.0 km", angle=angle
        )
