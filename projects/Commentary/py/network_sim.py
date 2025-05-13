import datetime as dt
import sys

sys.path.append("py/framework")

from intermagnet import Bfield
from mt_sites import PROFILES, TransferFunction
from network import Substation
from plots import StackPlots


def run_network_station_simulations_for_benchmark_event(
    site_name: str,
    stn: str,
    title: str = "",
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
    net = Substation(name=stn)
    gic_net = net.compute_J(ef.ex, ef.ey)

    sp = StackPlots(nrows=2, ncols=3)
    _, ax = sp.plot_stack_plots(
        bf.date,
        ef.ex / 1e3,
        dict(start=(0, 0), colspan=2, rowspan=1),
        text=f"E-field under {title}",
        ylabel="E-field (V/km)",
        color="r",
        lw=0.4,
        label=r"$E_x$",
        xlim=date_lim,
        ylim=[-2, 2],
        datetime=True,
        tag="(A)",
    )
    sp.plot_stack_plots(
        bf.date,
        ef.ey / 1e3,
        None,
        label=r"$E_y$",
        color="k",
        lw=0.4,
        xlim=date_lim,
        ax=ax,
        datetime=True,
        ylim=[-2, 2],
    )
    ax.legend(loc=1)
    _, ax = sp.plot_stack_plots(
        bf.date,
        gic_net / 1e3,
        dict(start=(1, 0), colspan=2, rowspan=1),
        text=rf"GIC through substation, {stn}",
        ylabel="GIC (A)",
        color="k",
        lw=0.9,
        xlim=date_lim,
        ylim=[-500, 500],
        xlabel="Time (UT)",
        datetime=True,
        tag="(B)",
    )

    theta, cor = net.compute_segmented_correlation(
        bf.bmag, ef.ex, ef.ey, normalize=site_name == "CaseA"
    )
    # sp = StackPlots(nrows=2, ncols=1, polar=True)
    sp.plot_dirctional_plots(
        theta,
        cor,
        dict(start=(0, 2), colspan=1, rowspan=1),
        title=f"{title} / Station: {stn}",
        text=r"r(GIC, $B_h$)",
        color="r",
        tag="(C)",
    )
    theta, cor = net.compute_segmented_correlation(
        bf.dbh, ef.ex, ef.ey, normalize=site_name == "CaseB"
    )
    sp.plot_dirctional_plots(
        theta,
        cor,
        dict(start=(1, 2), colspan=1, rowspan=1),
        text=r"r(GIC, $\partial B_h$)",
        tag="(D)",
        color="b",
    )
    sp.save_fig(
        f"figures/Comp_Paper_Benchmark_Net_{site_name}_{stn}.png"
    )
    sp.close()
    return

def compile_stack_plots(
    site_name: str,
    stn_names: list[str],
    title: str = "",
    date_lim: list[dt.datetime] = [
        dt.datetime(1989, 3, 12, 12),
        dt.datetime(1989, 3, 15),
    ],
):
    site = getattr(PROFILES, site_name)
    bf = Bfield.create_benchmark_bfield()
    tf = TransferFunction(site)
    ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta)
    sp = StackPlots(nrows=len(stn_names)+1, ncols=2)
    _, ax = sp.plot_stack_plots(
        bf.date,
        ef.ex / 1e3,
        dict(start=(0, 0), colspan=2, rowspan=1),
        text=f"E-field under {title}",
        ylabel="E-field (V/km)",
        color="r",
        lw=0.4,
        label=r"$E_x$",
        xlim=date_lim,
        ylim=[-2, 2],
        datetime=True,
        tag="(A)",
    )
    sp.plot_stack_plots(
        bf.date,
        ef.ey / 1e3,
        None,
        label=r"$E_y$",
        color="k",
        lw=0.4,
        xlim=date_lim,
        ax=ax,
        datetime=True,
        ylim=[-2, 2],
    )
    ax.legend(loc=1)

    for d, stn, tag in zip(
        range(len(stn_names)),
        stn_names, 
        [
            "(B)", "(C)", "(D)", 
            "(E)", "(F)", "(G)"
        ]
    ):
        net = Substation(name=stn)
        gic_net = net.compute_J(ef.ex, ef.ey)
        _, ax = sp.plot_stack_plots(
            bf.date,
            gic_net / 1e3,
            dict(start=(d+1, 0), colspan=2, rowspan=1),
            text=rf"GIC through substation, {stn}",
            ylabel="GIC (A)",
            color="k",
            lw=0.9,
            xlim=date_lim,
            ylim=[-500, 500],
            xlabel="Time (UT)" if d==len(stn_names)-1 else "",
            datetime=True,
            tag=tag,
        )
    sp.save_fig(
        f"figures/Stack_Paper_Benchmark_Net_{site_name}.png"
    )
    sp.close()
    return

def dial_plots(
    site_name: str,
    stn_names: list[str],
    title: str = "",
    date_lim: list[dt.datetime] = [
        dt.datetime(1989, 3, 12, 12),
        dt.datetime(1989, 3, 15),
    ],
):
    site = getattr(PROFILES, site_name)
    bf = Bfield.create_benchmark_bfield()
    tf = TransferFunction(site)
    ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta)
    sp = StackPlots(nrows=4, ncols=3)

    for d, stn, tag in zip(
        range(len(stn_names)),
        stn_names, 
        [
            "A", "B", "C", 
            "D", "E", "F"
        ]
    ):
        net = Substation(name=stn)

        theta, cor = net.compute_segmented_correlation(
            bf.bmag, ef.ex, ef.ey, normalize=site_name == "CaseA"
        )
        sp.plot_dirctional_plots(
            theta,
            cor,
            dict(start=(0+2*int(d/3), d%3), colspan=1, rowspan=1),
            title=f"{title} / Station: {stn}" if d==0 else f"Station: {stn}",
            text=r"r(GIC, $B_h$)",
            color="r",
            tag=f"({tag}-1)",
        )
        theta, cor = net.compute_segmented_correlation(
            bf.dbh, ef.ex, ef.ey, normalize=site_name == "CaseB"
        )
        sp.plot_dirctional_plots(
            theta,
            cor,
            dict(start=(1+2*int(d/3), d%3), colspan=1, rowspan=1),
            text=r"r(GIC, $\partial B_h$)",
            tag=f"({tag}-2)",
            color="b",
        )
    sp.save_fig(
        f"figures/Dial_Paper_Benchmark_Net_{site_name}.png"
    )
    sp.close()
    return


if __name__ == "__main__":
    dial_plots(
        "CaseA",
        ["S2", "S3", "S4", "S5", "S6", "S8"],
        r"Case A, $\tau_1$=10.0 km",
    )
    dial_plots(
        "CaseB",
        ["S2", "S3", "S4", "S5", "S6", "S8"],
        r"Case B, $\tau_1$=100.0 km",
    )
    compile_stack_plots(
        "CaseA",
        ["S2", "S3", "S4", "S5", "S6", "S8"],
        r"Case A, $\tau_1$=10.0 km",
    )
    compile_stack_plots(
        "CaseB",
        ["S2", "S3", "S4", "S5", "S6", "S8"],
        r"Case B, $\tau_1$=100.0 km",
    )
    
    # for stn in ["S2", "S3", "S4", "S5", "S6", "S8"]:
    #     run_network_station_simulations_for_benchmark_event(
    #         "CaseA",
    #         stn,
    #         r"Case A, $\tau_1$=10.0 km",
    #     )
    #     run_network_station_simulations_for_benchmark_event(
    #         "CaseB",
    #         stn,
    #         r"Case A, $\tau_1$=100.0 km",
    #     )
