import datetime as dt
import sys

sys.path.append("py/framework")

import numpy as np
from intermagnet import Bfield
from mt_sites import PROFILES, TransferFunction
from network import Substation
from plots import StackPlots


def run_network_stations_simulations_for_benchmark_event(
    site_name: str,
    stn_names: list[str],
    title: str = "",
    date_lim: list[dt.datetime] = [
        dt.datetime(1989, 3, 12, 12),
        dt.datetime(1989, 3, 15),
    ],
):
    """
    Run the pipeline simulation for the benchmark event.
    (No polar plots in this function — no changes needed to correlation calls.)
    """
    site = getattr(PROFILES, site_name)
    bf = Bfield.create_benchmark_bfield()
    tf = TransferFunction(site)
    ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta)

    sp = StackPlots(nrows=7, ncols=1)
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
        tag="(a)",
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

    for j, stn in enumerate(stn_names):
        net = Substation(name=stn)
        gic_net = net.compute_J(ef.ex, ef.ey)
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
            datetime=True,
            tag=f"({chr(98 + j)})",
        )
    sp.save_fig(
        f"figures/Comp_Paper_Benchmark_Net_{site_name}.png"
    )
    sp.close()
    return

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
    Run the station simulation for the benchmark event.
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

    # --- CHANGED: bf.bmag -> bf.bx, bf.by; removed normalize ---
    theta, cor = net.compute_segmented_correlation(
        bf.bx, bf.by, ef.ex, ef.ey
    )
    sp.plot_dirctional_plots(
        theta,
        np.abs(cor),
        dict(start=(0, 2), colspan=1, rowspan=1),
        title=f"{title} / Station: {stn}",
        text=r"r(GIC, $B_\theta$)",
        color="r",
        tag="(C)",
    )
    # --- CHANGED: bf.dbh -> bf.dbx, bf.dby; removed normalize ---
    theta, cor = net.compute_segmented_correlation(
        bf.dbx, bf.dby, ef.ex, ef.ey
    )
    sp.plot_dirctional_plots(
        theta,
        np.abs(cor),
        dict(start=(1, 2), colspan=1, rowspan=1),
        text=r"r(GIC, $\frac{\partial B_\theta}{\partial t}$)",
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
    """No correlation calls here — no changes needed."""
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
        ["(B)", "(C)", "(D)", "(E)", "(F)", "(G)"]
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
    stn_names: list[str],
):
    siteA, siteB = getattr(PROFILES, "CaseA"),getattr(PROFILES, "CaseB")
    bf = Bfield.create_benchmark_bfield()
    tfA, tfB = TransferFunction(siteA), TransferFunction(siteB)
    efA, efB = (
        tfA.compute_Efield(bf.bx, bf.by, bf.del_ta),
        tfB.compute_Efield(bf.bx, bf.by, bf.del_ta)
    )
    sp = StackPlots(nrows=6, ncols=2)

    for d, stn, tag in zip(
        range(len(stn_names)),
        stn_names, 
        ["A", "B", "C", "D", "E", "F"]
    ):
        net = Substation(name=stn)

        # --- CHANGED: bf.bmag -> bf.bx, bf.by; removed normalize ---
        thetaAb, corAb = net.compute_segmented_correlation(
            bf.bx, bf.by, efA.ex, efA.ey
        )
        # --- CHANGED: bf.dbh -> bf.dbx, bf.dby; removed normalize ---
        thetaAdb, corAdb = net.compute_segmented_correlation(
            bf.dbx, bf.dby, efA.ex, efA.ey
        )
        _, ax = sp.plot_dirctional_plots(
            thetaAb,
            np.abs(corAb),  # --- CHANGED: added np.abs() ---
            dict(start=(d, 0), colspan=1, rowspan=1),
            title=r"Case A, $\tau_1$=10.0 km" if d == 0 else None,
            text=r"r(GIC, $B_\theta$)" if d==0 else None,
            color="r",
            tag=f"({tag}-1)",
        )
        ax.text(-0.4, 0.5, f"Station: {stn}", ha="left", va="center", transform=ax.transAxes, rotation=90)
        sp.plot_dirctional_plots(
            thetaAdb,
            np.abs(corAdb),  # --- CHANGED: added np.abs() ---
            dict(start=(d, 0), colspan=1, rowspan=1),
            tag=f"",
            color="b",
            ax=ax,
        )

        # --- CHANGED: bf.bmag -> bf.bx, bf.by; removed normalize ---
        thetaBb, corBb = net.compute_segmented_correlation(
            bf.bx, bf.by, efB.ex, efB.ey
        )
        # --- CHANGED: bf.dbh -> bf.dbx, bf.dby; removed normalize ---
        thetaBdb, corBdb = net.compute_segmented_correlation(
            bf.dbx, bf.dby, efB.ex, efB.ey
        )
        _, ax = sp.plot_dirctional_plots(
            thetaBb,
            np.abs(corBb),  # --- CHANGED: added np.abs() ---
            dict(start=(d, 1), colspan=1, rowspan=1),
            title=r"Case B, $\tau_1$=100.0 km" if d == 0 else None,
            color="r",
            tag=f"({tag}-2)",
        )
        sp.plot_dirctional_plots(
            thetaBdb,
            np.abs(corBdb),  # --- CHANGED: added np.abs() ---
            dict(start=(d, 1), colspan=1, rowspan=1),
            text=r"r(GIC, $\frac{\partial B_\theta}{\partial t}$)" if d == 0 else None,
            tag=f"",
            color="b",
            ax=ax,
        )
    sp.fig.subplots_adjust(wspace=0.1)
    sp.save_fig(
        f"figures/Dial_Paper_Benchmark_Net.png"
    )
    sp.close()
    return


if __name__ == "__main__":
    dial_plots(
        ["S2", "S3", "S4", "S5", "S6", "S8"],
    )
    
    run_network_stations_simulations_for_benchmark_event(
        "CaseA",
        ["S2", "S3", "S4", "S5", "S6", "S8"],
        r"Case A, $\tau_1$=10.0 km",
    )
