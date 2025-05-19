import datetime as dt

import numpy as np
import pandas as pd
from scipy.signal import detrend

convertable_parameters = [
    dict(
        base_folder="10-11Feb1958/",
        file_details=[
            dict(
                component="Voltage",
                file="TAT1Volt.csv",
                y1=2000,
                y2=500,
                ybase=0,
                yscale=1,
                x1=4,
                x2=1,
                x1_date=dt.datetime(1958, 2, 11, 4),
                x2_date=dt.datetime(1958, 2, 11, 1),
                xlim=[dt.datetime(1958, 2, 11, 1), dt.datetime(1958, 2, 11, 4)],
                ylabel=r"TAT-1 Voltage [V]",
                text="Scaled Voltage from TAT1, 11 Feb 1958",
                ylim=[-3000, 3000],
                scale_details=dict(
                    xs=[1, 0.889, 0.811, 0.785],
                    ys=[0, 500, 1000, 1500],
                ),
            ),
        ],
    )
]


def rescale_data(Y, ybase, yscale):
    d = detrend(Y, type="constant")
    d = ybase + (d * yscale)
    return d


def rescale_time(T, x1, x2, t1, t2):
    dtau = (t2 - t1).total_seconds()
    tscale = [dtau * (t - x1) / (x2 - x1) for t in T]
    t = [t1 + dt.timedelta(seconds=s) for s in tscale]
    return t, tscale


def draw_dataset(scaled, param, ylabel, fname, xlim, text, fitted=None, ylim=None):
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.dates import DateFormatter
    import scienceplots
    plt.style.use(["science", "ieee"])
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [
        "Tahoma",
        "DejaVu Sans",
        "Lucida Grande",
        "Verdana",
    ]
    plt.rcParams.update({
        "text.usetex": False,
    })

    fig = plt.figure(dpi=300, figsize=(6, 3))
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(DateFormatter(r"$%H^{%M}$"))
    hours = mdates.HourLocator(byhour=range(0, 24, 1))
    ax.xaxis.set_major_locator(hours)
    ax.plot(scaled.Time, scaled[param], "ro", ms=0.5, ls="None")
    if fitted:
        ax.plot(fitted.Time, fitted[param], "b", lw=0.6, ls="-")
    ax.set_xlabel(r"Time [UT]")
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ylim = (
        ylim
        if ylim
        else [
            int(np.min(scaled[param]) / 1000) * 1e3,
            int(np.max(scaled[param]) / 1000) * 1e3 + 1000,
        ]
    )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.text(0.1, 0.9, text, ha="left", va="center", transform=ax.transAxes)
    fig.savefig(fname, bbox_inches="tight")
    return

def recalibrate_X(scale_details, Time, Volt):
    xs, ys = (
        np.array(scale_details["xs"]), 
        np.array(scale_details["ys"])
    )
    xds, yds = (
        np.abs(xs-xs.max()),
        np.abs(ys-ys.max())
    ) 
    coef = np.polyfit(yds, xds, 2)
    yds_new = np.abs(np.array(Volt) - ys.max())
    xds_new = np.polyval(coef, yds_new)
    raw = pd.DataFrame()
    raw["Time"], raw["Voltage"] = Time + xds_new, Volt
    return raw

def rebound_Voltage(o):
    o.Voltage *= 2.7 
    return o


for cp in convertable_parameters:
    for fd in cp["file_details"]:
        o = pd.DataFrame()
        raw = pd.read_csv(cp["base_folder"] + f"{fd['component']}/" + fd["file"])
        raw = recalibrate_X(fd["scale_details"], raw.Time, raw.Voltage)
        o[fd["component"]] = rescale_data(
            np.array(raw[fd["component"]]), fd["ybase"], fd["yscale"]
        )
        o["Time"], o["dT"] = rescale_time(
            np.array(raw["Time"]), fd["x1"], fd["x2"], fd["x1_date"], fd["x2_date"]
        )
        o = rebound_Voltage(o)
        o.to_csv(
            cp["base_folder"]
            + f"{fd['component']}/"
            + fd["file"].replace(".csv", "-rescale.csv"),
            index=False,
            header=True,
        )
        o = o.sort_values(by="Time")
        draw_dataset(
            o,
            param=fd["component"],
            ylabel=fd["ylabel"],
            fname=cp["base_folder"]
            + f"{fd['component']}/"
            + fd["file"].replace(".csv", "-scaled.png"),
            xlim=fd["xlim"],
            text=fd["text"],
            ylim=fd["ylim"],
        )
        o.Time.iloc[0] = o.Time.iloc[0].to_pydatetime().replace(second=0, microsecond=0)
        h = o.set_index("Time")
        h = h[~h.index.duplicated()]
        new_index_20s = pd.date_range(h.index.min(), h.index.max(), freq="1s")
        tmp_index_20s = h.index.union(new_index_20s)
        df = h.reindex(tmp_index_20s).interpolate("cubic").reindex(new_index_20s)
        df.index.name = "Time"
        df = df.reset_index()
        print(df.head())
        draw_dataset(
            df,
            param=fd["component"],
            ylabel=fd["ylabel"],
            fname=cp["base_folder"]
            + f"{fd['component']}/"
            + fd["file"].replace(".csv", "-scaled-HR.png"),
            xlim=fd["xlim"],
            text=fd["text"],
            ylim=fd["ylim"],
        )
        df.to_csv(
            cp["base_folder"]
            + f"{fd['component']}/"
            + fd["file"].replace(".csv", "-rescale-HR.csv"),
            index=False,
            header=True,
        )
