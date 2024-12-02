import numpy as np
import pandas as pd
import datetime as dt

from scipy.signal import detrend

convertable_parameters = [
    dict(
        base_folder = "10-11Feb1958/",
        file_details = [
            dict(
                component="H",
                file="H[Eskdalemuir].csv",
                y1=80, y2=95, ybase=15894, yscale=16.6,
                x1=60, x2=275, 
                x1_date=dt.datetime(1958,2,10,12),
                x2_date=dt.datetime(1958,2,11,2),
                xlim=[
                    dt.datetime(1958,2,10,18),
                    dt.datetime(1958,2,11,6)
                ],
                ylabel=r"H-Component [$\gamma$]",
                text="Scaled H-Component from Eskdalemuir, 10-11 Feb 1958",
                ylim=None,
            ),
            dict(
                component="D",
                file="D[Eskdalemuir].csv",
                y1=75, y2=92, ybase=9.8675, yscale=0.016667,
                x1=38, x2=237, 
                x1_date=dt.datetime(1958,2,10,12),
                x2_date=dt.datetime(1958,2,11),
                xlim=[
                    dt.datetime(1958,2,10,18),
                    dt.datetime(1958,2,11,6)
                ],
                ylabel=r"D-Component [$\circ$]",
                text="Scaled D-Component from Eskdalemuir, 10-11 Feb 1958",
                ylim=[8, 13]
            ),
            dict(
                component="Z",
                file="Z[Eskdalemuir].csv",
                y1=103, y2=110, ybase=44765, yscale=8.65,
                x1=58, x2=245, 
                x1_date=dt.datetime(1958,2,10,12),
                x2_date=dt.datetime(1958,2,11),
                xlim=[
                    dt.datetime(1958,2,10,18),
                    dt.datetime(1958,2,11,6)
                ],
                ylabel=r"Z-Component [$\gamma$]",
                text="Scaled Z-Component from Eskdalemuir, 10-11 Feb 1958",
                ylim=[43000, 46000]
            )
        ]
    )
]

def rescale_data(Y, ybase, yscale):
    d = detrend(Y, type="constant")
    d = ybase + (d * yscale)
    return d

def rescale_time(T, x1, x2, t1, t2):
    dtau = (t2-t1).total_seconds()
    tscale = [dtau*(t-x1)/(x2-x1) for t in T] 
    t = [t1 + dt.timedelta(seconds=s) for s in tscale]
    return t, tscale

def draw_dataset(scaled, param, ylabel, fname, xlim, text, fitted=None, ylim=None):
    import matplotlib.pyplot as plt
    from matplotlib.dates import DateFormatter
    import matplotlib.dates as mdates
    plt.style.use(["science", "ieee"])
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                    "Lucida Grande", "Verdana"]
    
    fig = plt.figure(dpi=300, figsize=(6, 3))
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(DateFormatter(r"%H^{%M}"))
    hours = mdates.HourLocator(byhour=range(0, 24, 1))
    ax.xaxis.set_major_locator(hours)
    ax.plot(scaled.Time, scaled[param], "ro", ms=0.5, ls="None")
    if fitted:
        ax.plot(fitted.Time, fitted[param], "b", lw=0.6, ls="-")
    ax.set_xlabel(r"Time [UT]")
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ylim = ylim if ylim else [
        int(np.min(scaled[param])/1000)*1e3,
        int(np.max(scaled[param])/1000)*1e3 + 1000,
    ]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.text(0.1, 0.9, text, ha="left", va="center", transform=ax.transAxes)
    fig.savefig(fname, bbox_inches="tight")
    return

for cp in convertable_parameters:
    for fd in cp["file_details"]:
        o = pd.DataFrame()
        raw = pd.read_csv(cp["base_folder"]+f"{fd['component']}/"+fd["file"])
        o[fd["component"]] = rescale_data(
            np.array(raw[fd["component"]]),
            fd["ybase"], fd["yscale"]
        )
        o["Time"], o["dT"] = rescale_time(
            np.array(raw["Time"]), fd["x1"], fd["x2"],
            fd["x1_date"], fd["x2_date"]
        )
        o.to_csv(
            cp["base_folder"] +f"{fd['component']}/"+ fd["file"].replace(".csv", "-rescale.csv"), 
            index=False, header=True
        )
        o = o.sort_values(by="Time")
        draw_dataset(
            o,
            param=fd["component"], ylabel=fd["ylabel"],
            fname=cp["base_folder"] +f"{fd['component']}/"+ fd["file"].replace(".csv", "-scaled.png"),
            xlim=fd["xlim"], text=fd["text"], ylim=fd["ylim"]
        )
        o.Time.iloc[0] = o.Time.iloc[0].to_pydatetime().replace(second=0,microsecond=0)
        h = o.set_index("Time")
        h = h[~h.index.duplicated()]
        new_index_20s = pd.date_range(h.index.min(), h.index.max(), freq='20s')
        tmp_index_20s = h.index.union(new_index_20s)
        df = h.reindex(tmp_index_20s).interpolate("index").reindex(new_index_20s)
        df.index.name = "Time"
        df = df.reset_index()
        print(df.head())
        draw_dataset(
            df,
            param=fd["component"], ylabel=fd["ylabel"],
            fname=cp["base_folder"] +f"{fd['component']}/"+ fd["file"].replace(".csv", "-scaled-HR.png"),
            xlim=fd["xlim"], text=fd["text"], ylim=fd["ylim"]
        )
        df.to_csv(
            cp["base_folder"] +f"{fd['component']}/"+ fd["file"].replace(".csv", "-rescale-HR.csv"), 
            index=False, header=True
        )