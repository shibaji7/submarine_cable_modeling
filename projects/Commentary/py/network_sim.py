import datetime as dt
import sys

sys.path.append("py/framework")

from intermagnet import Bfield
from mt_sites import PROFILES, TransferFunction
from network import Substation
from plots import StackPlots

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(["science", "ieee"])
plt.rcParams.update(
    {
        "text.usetex": True,
    }
)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]
size = 18
mpl.rcParams.update(
    {"xtick.labelsize": size, "ytick.labelsize": size, "font.size": size}
)


def run_network_station_simulations_for_benchmark_event(
    site_name: str,
    stn: str,
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

    # --- CHANGED: pass bx, by instead of bf.bmag ---
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

    # --- CHANGED: pass dbx, dby instead of bf.dbh ---
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


def compile_double_stack_plots(
    site_names: list[str],
    stn_names: list[str],
    titles: list[str] = ["", ""],
    date_lim: list[dt.datetime] = [
        dt.datetime(1989, 3, 12, 12),
        dt.datetime(1989, 3, 15),
    ],
):
    bf = Bfield.create_benchmark_bfield()
    sp = StackPlots(nrows=len(stn_names)+1, ncols=4)
    dpx = 0
    for site_name, title in zip(site_names, titles):
        site = getattr(PROFILES, site_name)
        tf = TransferFunction(site)
        ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta)
        _, ax = sp.plot_stack_plots(
            bf.date,
            ef.ex / 1e3,
            dict(start=(0, 2*dpx), colspan=2, rowspan=1),
            text=f"({dpx+1}) E-field under {title}",
            ylabel="E-field (V/km)" if dpx == 0 else "",
            color="r",
            lw=0.4,
            label=r"$E_x$",
            xlim=date_lim,
            ylim=[-2, 2],
            datetime=True,
            tag="(a-" + str(dpx+1) + ")",
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
        ax.set_xticklabels([])
        ax.set_xticklabels([], minor=True)
        if dpx > 0:
            ax.set_yticklabels([])

        for d, stn, tag in zip(
            range(len(stn_names)),
            stn_names, 
            ["b", "c", "d", "e", "f", "g"]
        ):
            net = Substation(name=stn)
            gic_net = net.compute_J(ef.ex, ef.ey)
            _, ax = sp.plot_stack_plots(
                bf.date,
                gic_net / 1e3,
                dict(start=(d+1, 2*dpx), colspan=2, rowspan=1),
                text="",
                ylabel="GIC (A)" if dpx == 0 else "",
                color="k",
                lw=0.9,
                xlim=date_lim,
                ylim=[-500, 500],
                xlabel="Time (UT)" if d==len(stn_names)-1 else "",
                datetime=True,
                tag=f"({tag}-{dpx+1}) GIC through substation, {stn}",
            )
            if d != len(stn_names) - 1:
                ax.set_xticklabels([])
                ax.set_xticklabels([], minor=True)
            if dpx > 0:
                ax.set_yticklabels([])
        dpx += 1
    sp.fig.subplots_adjust(hspace=0.0, wspace=0.2)
    sp.save_fig(
        f"figures/Stack_Paper_Benchmark_Net.png"
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
    ax.set_xticklabels([])
    ax.set_xticklabels([], minor=True)

    for d, stn, tag in zip(
        range(len(stn_names)),
        stn_names, 
        ["(b)", "(c)", "(d)", "(e)", "(f)", "(g)"]
    ):
        net = Substation(name=stn)
        gic_net = net.compute_J(ef.ex, ef.ey)
        _, ax = sp.plot_stack_plots(
            bf.date,
            gic_net / 1e3,
            dict(start=(d+1, 0), colspan=2, rowspan=1),
            text="",
            ylabel="GIC (A)",
            color="k",
            lw=0.9,
            xlim=date_lim,
            ylim=[-500, 500],
            xlabel="Time (UT)" if d==len(stn_names)-1 else "",
            datetime=True,
            tag=tag + rf" GIC through substation, {stn}",
        )
        if d != len(stn_names) - 1:
            ax.set_xticklabels([])
            ax.set_xticklabels([], minor=True)
    sp.fig.subplots_adjust(hspace=0.0, wspace=0.0)
    sp.save_fig(
        f"figures/Stack_Paper_Benchmark_Net_{site_name}.png"
    )
    sp.close()
    return


def dial_double_plots(
    site_names: list[str],
    stn_names: list[str],
    titles: list[str] = ["", ""],
    date_lim: list[dt.datetime] = [
        dt.datetime(1989, 3, 12, 12),
        dt.datetime(1989, 3, 15),
    ],
):
    bf = Bfield.create_benchmark_bfield()
    sp = StackPlots(nrows=6, ncols=2)
    dpx = 0
    for site_name, title in zip(site_names, titles):
        site = getattr(PROFILES, site_name)
        tf = TransferFunction(site)
        ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta)
        
        for d, stn, tag in zip(
            range(len(stn_names)),
            stn_names, 
            ["a", "b", "c", "d", "e", "f"]
        ): 
            print(f"Processing {site_name} - {stn}...")
            net = Substation(name=stn)

            # --- CHANGED: pass bx, by instead of bf.bmag ---
            theta, cor = net.compute_segmented_correlation(
                bf.bx, bf.by, ef.ex, ef.ey
            )
            title_txt = f"Substation: {stn}" if dpx == 0 else ""
            _, ax = sp.plot_dirctional_plots(
                theta,
                np.abs(cor),
                dict(start=(d, dpx), colspan=1, rowspan=1),
                title=f"({dpx+1}) {title}" if d==0 else "",
                text=r"r(GIC, $B_\theta$)" if d==2 and dpx==0 else "",
                color="r",
                tag=f"({tag}-{dpx+1})",
                text_location=(1.2, 1.2)
            )
            ax.text(-0.23, 0.5, title_txt, transform=ax.transAxes, ha="center", va="center", rotation=90)

            # --- CHANGED: pass dbx, dby instead of bf.dbh ---
            theta, cor = net.compute_segmented_correlation(
                bf.dbx, bf.dby, ef.ex, ef.ey
            )
            sp.plot_dirctional_plots(
                theta,
                np.abs(cor),
                None,
                text=r"r(GIC, $\frac{\partial B_\theta}{\partial t}$)" if d==4 and dpx==0 else "",
                color="b",
                ax=ax,
                text_location=(1.2, 1.2)
            )
        dpx += 1
    sp.fig.subplots_adjust(hspace=0.3, wspace=0.2)
    sp.save_fig(
        f"figures/Dial_Paper_Benchmark_Net.png"
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
    sp = StackPlots(nrows=2, ncols=3)

    for d, stn, tag in zip(
        range(len(stn_names)),
        stn_names, 
        ["a", "b", "c", "d", "e", "f"]
    ):
        net = Substation(name=stn)

        # --- CHANGED: pass bx, by instead of bf.bmag ---
        theta, cor = net.compute_segmented_correlation(
            bf.bx, bf.by, ef.ex, ef.ey
        )
        _, ax = sp.plot_dirctional_plots(
            theta,
            np.abs(cor),
            dict(start=(0+int(d/3), d%3), colspan=1, rowspan=1),
            title=f"{title} / Substation: {stn}" if d==0 else f"Substation: {stn}",
            text=r"r(GIC, $B_\theta$)" if d==4 else "",
            color="r",
            tag=f"({tag})",
            text_location=(-.65, 1.2)
        )

        # --- CHANGED: pass dbx, dby instead of bf.dbh ---
        theta, cor = net.compute_segmented_correlation(
            bf.dbx, bf.dby, ef.ex, ef.ey
        )
        sp.plot_dirctional_plots(
            theta,
            np.abs(cor),
            None,
            text=r"r(GIC, $\frac{\partial B_\theta}{\partial t}$)" if d==4 else "",
            color="b",
            ax=ax,
            text_location=(1.2, 1.2)
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
        r"Case A",
    )
    dial_plots(
        "CaseB",
        ["S2", "S3", "S4", "S5", "S6", "S8"],
        r"Case B",
    )
    compile_stack_plots(
        "CaseA",
        ["S2", "S3", "S4", "S5", "S6", "S8"],
        "Case A",
    )
    compile_stack_plots(
        "CaseB",
        ["S2", "S3", "S4", "S5", "S6", "S8"],
        "Case B",
    )

    compile_double_stack_plots(
        ["CaseA", "CaseB"],
        ["S2", "S3", "S4", "S5", "S6", "S8"],
        titles=["Case A", "Case B"],
    )
    dial_double_plots(
        ["CaseA", "CaseB"],
        ["S2", "S3", "S4", "S5", "S6", "S8"],
        ["Case A", "Case B"],
    )
