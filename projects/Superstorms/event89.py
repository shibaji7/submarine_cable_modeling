import os
import pandas as pd
import datetime as dt
import numpy as np
from loguru import logger

os.environ["OMNIDATA_PATH"] = "/home/shibaji/omni/"

def _load_omni_(dates, res=1):
    import pyomnidata
    logger.info(f"OMNIDATA_PATH: {os.environ['OMNIDATA_PATH']}")
    pyomnidata.UpdateLocalData()
    omni = pd.DataFrame(
        pyomnidata.GetOMNI(dates[0].year,Res=res)
    )
    omni["time"] = omni.apply(
        lambda r: (
            dt.datetime(
                int(str(r.Date)[:4]), 
                int(str(r.Date)[4:6]),
                int(str(r.Date)[6:].replace(".0","")) 
            ) 
            + dt.timedelta(hours=r.ut)
        ), 
        axis=1
    )
    omni = omni[
        (omni.time>=dates[0])
        & (omni.time<=dates[1])
    ]
    return omni

def create_stack_plots(dates, size=15):
    omni = _load_omni_(dates)
    
    import matplotlib.dates as mdates
    from matplotlib.dates import DateFormatter

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.style.use(["science", "ieee"])
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                    "Lucida Grande", "Verdana"]
    mpl.rcParams.update(
        {"xtick.labelsize": size, "ytick.labelsize": size, "font.size": size}
    )
    
    fig, axes = plt.subplots(nrows=2, ncols=1, dpi=300, figsize=(6, 4), sharex=True)

    # ax = axes[0]
    # ax.xaxis.set_major_formatter(DateFormatter(r"%H UT"))
    # ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
    # ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 1)))
    # ax.plot(
    #     omni.time, omni.BzGSM, 
    #     color="r", ls="-", lw=0.8, label=r"$B_z$"
    # )
    # ax.plot(
    #     omni.time, omni.ByGSM, 
    #     color="g", ls="-", lw=0.8, label=r"$B_y$"
    # )
    # ax.plot(
    #     omni.time, omni.BxGSE, 
    #     color="k", ls="-", lw=0.8, label=r"$B_x$"
    # )
    # ax.set_ylim(-100, 100)
    # ax.set_xlim(dates)
    # ax.legend(loc=1)
    # ax.set_ylabel("IMF, nT")

    # ax = axes[1]
    # ax.xaxis.set_major_formatter(DateFormatter(r"%H UT"))
    # ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
    # ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 1)))
    # ax.plot(
    #     omni.time, omni.ProtonDensity, 
    #     color="k", ls="-", lw=0.8, label=r"n, /cc"
    # )
    # ax.plot(
    #     omni.time, omni.FlowPressure, 
    #     color="b", ls="-", lw=0.8, label=r"$P_{dyn}$"
    # )
    # ax.plot(
    #     omni.time, omni.FlowSpeed/10, 
    #     color="m", ls="-", lw=0.8, label=r"$V\times 10$, km/s"
    # )
    # ax.legend(loc=1)
    # ax.set_ylabel("SW Params")
    # ax.set_ylim(0, 120)
    # ax.set_xlim(dates)


    ax = axes[0]
    ax.xaxis.set_major_formatter(DateFormatter(r"%H UT"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 1)))
    ax.text(0.05, 1.05, "Date: 13-14 March 1989", ha="left", va="bottom", transform=ax.transAxes)
    ax.plot(omni.time, omni.SymH, "k", lw=0.9)
    ax.set_ylabel("SYM-H (nT)")
    ax.set_ylim(-800, 200)
    ax = ax.twinx()
    ax.plot(omni.time, omni.AsyH, "b", lw=0.8)
    ax.set_ylim(0, 1200)
    ax.set_ylabel("ASY-H (nT)", fontdict={"color":"b"})
    ax.set_xlim(dates)

    ax = axes[1]
    ax.plot(omni.time, omni.AU, "r-", lw=0.8, label="AU")
    ax.plot(omni.time, omni.AL, "b-", lw=0.8, label="AL")
    ax.plot(omni.time, omni.AE, "k-", lw=0.8, label="AE")
    ax.legend(loc=1)
    ax.set_ylim(-6000, 6000)
    # ax.set_yticks([-4, -2, 0, 2, 4])
    # ax.set_yticklabels([r"-$10^{4}$", r"-$10^{2}$", 0, r"$10^{2}$", r"$10^{4}$"])
    ax.set_ylabel("AE/AL/AU (nT)")
    ax.set_xlim(dates)

    for ax in axes:
        ax.axvline(dt.datetime(1989,3,13,1,30), ls="--", lw=1.5, color="r")
        ax.axvline(dt.datetime(1989,3,13,11,10), ls="--", lw=1.5, color="k")
        ax.axvline(dt.datetime(1989,3,13,21,45), ls="--", lw=1.5, color="g")
        ax.axvline(dt.datetime(1989,3,14,1,30), ls="--", lw=1.5, color="b")
    print(omni.columns)
    ax.set_xlabel("Time [UT]")
    fig.savefig("figures/Event.png", bbox_inches="tight")
    return

if __name__ == "__main__":
    dates = [dt.datetime(1989,3,13), dt.datetime(1989,3,14,12)]
    create_stack_plots(dates)