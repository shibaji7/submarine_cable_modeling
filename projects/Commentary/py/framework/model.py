import datetime as dt
import sys

import numpy as np

sys.path.append("framework")

from intermagnet import Bfield
from mt_sites import PROFILES, TransferFunction
from network import Substation
from plots import StackPlots


def run_analysis(site, title):
    delt = 24 * 6 * 60
    dxlim = [dt.datetime(1989, 3, 12, 12), dt.datetime(1989, 3, 15)]
    bf = Bfield.create_benchmark_bfield()
    tf = TransferFunction(site)
    ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta, p=0.0)

    B, label, color, tlab = bf.by, r"$B_y$", "k", "$B_y$"
    (s2, s3, s4, s5, s6, s8) = (
        Substation(name="S2"),
        Substation(name="S3"),
        Substation(name="S4"),
        Substation(name="S5"),
        Substation(name="S6"),
        Substation(name="S8"),
    )

    sp = StackPlots(nrows=4, ncols=1, datetime=True)
    _, ax = sp.plot_stack_plots(
        bf.date,
        B,
        text=title,
        ylabel="B-field (nT)",
        color=color,
        label=label,
        lw=0.4,
        xlim=dxlim,
        ylim=[-4000, 4000],
    )
    ax.legend(loc=1)

    _, ax = sp.plot_stack_plots(
        bf.date,
        s2.compute_J(ef.ex / 1e3, ef.ey / 1e3),
        text="",
        ylabel="GIC at S2 (A)",
        color="b",
        lw=0.4,
        xlim=dxlim,
        ylim=[-400, 400],
        ylabel_color="b",
    )
    ax.text(
        0.05,
        0.95,
        rf"r(GIC@S2, {tlab})=%.3f"
        % np.round(
            np.corrcoef(
                s2.compute_J(ef.ex / 1e3, ef.ey / 1e3)[delt:-delt], B[delt:-delt]
            )[0, 1],
            3,
        ),
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(color="b", size=12),
    )
    sp.plot_stack_plots(
        bf.date,
        s3.compute_J(ef.ex / 1e3, ef.ey / 1e3),
        text="",
        ylabel="GIC at S3 (A)",
        color="k",
        ax=ax.twinx(),
        lw=0.4,
        xlim=dxlim,
        ylim=[-400, 400],
    )
    ax.text(
        0.95,
        0.95,
        rf"r(GIC@S3, {tlab})=%.3f"
        % np.round(
            np.corrcoef(
                s3.compute_J(ef.ex / 1e3, ef.ey / 1e3)[delt:-delt], B[delt:-delt]
            )[0, 1],
            3,
        ),
        ha="right",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size=12),
    )
    _, ax = sp.plot_stack_plots(
        bf.date,
        s4.compute_J(ef.ex / 1e3, ef.ey / 1e3),
        text="",
        ylabel="GIC at S4 (A)",
        color="b",
        lw=0.4,
        xlim=dxlim,
        ylim=[-400, 400],
        ylabel_color="b",
    )
    ax.text(
        0.05,
        0.95,
        rf"r(GIC@S4, {tlab})=%.3f"
        % np.round(
            np.corrcoef(
                s4.compute_J(ef.ex / 1e3, ef.ey / 1e3)[delt:-delt], B[delt:-delt]
            )[0, 1],
            3,
        ),
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(color="b", size=12),
    )
    sp.plot_stack_plots(
        bf.date,
        s5.compute_J(ef.ex / 1e3, ef.ey / 1e3),
        text="",
        ylabel="GIC at S5 (A)",
        color="k",
        ax=ax.twinx(),
        lw=0.4,
        xlim=dxlim,
        ylim=[-400, 400],
    )
    ax.text(
        0.95,
        0.95,
        rf"r(GIC@S5, {tlab})=%.3f"
        % np.round(
            np.corrcoef(
                s5.compute_J(ef.ex / 1e3, ef.ey / 1e3)[delt:-delt], B[delt:-delt]
            )[0, 1],
            3,
        ),
        ha="right",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size=12),
    )
    _, ax = sp.plot_stack_plots(
        bf.date,
        s6.compute_J(ef.ex / 1e3, ef.ey / 1e3),
        text="",
        ylabel="GIC at S6 (A)",
        color="b",
        lw=0.4,
        xlim=dxlim,
        ylim=[-400, 400],
        ylabel_color="b",
    )
    ax.text(
        0.05,
        0.95,
        rf"r(GIC@S4, {tlab})=%.3f"
        % np.round(
            np.corrcoef(
                s6.compute_J(ef.ex / 1e3, ef.ey / 1e3)[delt:-delt], B[delt:-delt]
            )[0, 1],
            3,
        ),
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(color="b", size=12),
    )
    sp.plot_stack_plots(
        bf.date,
        s8.compute_J(ef.ex / 1e3, ef.ey / 1e3),
        text="",
        ylabel="GIC at S8 (A)",
        color="k",
        ax=ax.twinx(),
        lw=0.4,
        xlim=dxlim,
        ylim=[-400, 400],
    )
    ax.text(
        0.95,
        0.95,
        rf"r(GIC@S8, {tlab})=%.3f"
        % np.round(
            np.corrcoef(
                s8.compute_J(ef.ex / 1e3, ef.ey / 1e3)[delt:-delt], B[delt:-delt]
            )[0, 1],
            3,
        ),
        ha="right",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size=12),
    )
    sp.save_fig("../../figures/Benchmark_Network_Stack.png")
    sp.close()
    return


def run_r_analysis(site, stn):
    delt = 10 * 6 * 60
    dxlim = [dt.datetime(1989, 3, 12, 12), dt.datetime(1989, 3, 15)]
    bf = Bfield.create_benchmark_bfield()
    tf = TransferFunction(site)
    ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta, p=0.0)

    B, label, color, tlab = bf.bx, r"$B_x$", "r", "$B_x$"
    s = Substation(name=stn)
    gic = s.compute_J(ef.ex / 1e3, ef.ey / 1e3)
    sp = StackPlots(nrows=1, ncols=1, datetime=False, figsize=(4, 4))
    ax = sp.axes[0]
    r = np.corrcoef(gic[delt:-delt], bf.bx[delt:-delt])[0, 1]
    ax.plot(
        gic[delt:-delt],
        bf.bx[delt:-delt],
        "r.",
        label=rf"r(GIC[@{stn}], $B_x$)=%.3f" % r,
    )
    r = np.corrcoef(gic[delt:-delt], bf.by[delt:-delt])[0, 1]
    ax.plot(
        gic[delt:-delt],
        bf.by[delt:-delt],
        "k.",
        label=rf"r(GIC[@{stn}], $B_y$)=%.3f" % r,
    )
    ax.set_xlabel("GIC, A")
    ax.set_ylabel("B, nT")
    ax.legend(loc=1)
    ax.set_ylim(-4000, 4000)
    sp.save_fig(f"../../figures/{stn}.png")
    sp.close()
    return


def run_network_station_simulations_for_benchmark_event(site, stn, title):
    """
    Run the station simulation for the benchmark event.
    """
    dxlim = [dt.datetime(1989, 3, 12, 12), dt.datetime(1989, 3, 15)]
    bf = Bfield.create_benchmark_bfield()
    tf = TransferFunction(site)
    ef = tf.compute_Efield(bf.bx, bf.by, bf.del_ta, p=0)
    station = Substation(name=stn)

    theta, cor = station.compute_segmented_correlation(bf.bmag, ef.ex, ef.ex, True)
    # theta, gic = station.compute_static_gic()
    sp = StackPlots(nrows=1, ncols=1, polar=True)
    sp.plot_dirctional_plots(
        theta,
        cor,
        title=f"{title} / Station: {stn}",
        text=r"r(GIC, B)",
        color="r",
        cable_angle=None,
    )
    # sp.plot_dirctional_plots(
    #     theta,
    #     np.abs(gic),
    #     title=f"{title} / Station: {stn}",
    #     text=r"GIC",
    #     color="r",
    #     rlims=[0, 500],
    #     rticks=[0, 200, 500],
    #     cable_angle=None,
    # )
    # theta, cor = pipeline.compute_segmented_correlation(bf.dbx, bf.dby, 0.5)
    # sp.plot_dirctional_plots(
    #     theta,
    #     cor,
    #     text=r"r(GIC, $\partial B_h$)",
    #     color="b",
    #     cable_angle=angle,
    # )
    sp.save_fig(f"../../figures/Benchmark_Network_Dirc_{stn}.png")
    sp.close()
    return


if __name__ == "__main__":
    run_network_station_simulations_for_benchmark_event(
        PROFILES.CaseA, "S2", r"Case A, $\tau_1$=10.0 km"
    )
    # run_analysis(PROFILES.CaseA, r"Case A, $\tau_1$=100.0 km")
    # run_r_analysis(PROFILES.CaseA, "S8")
