import pandas as pd
from scubas.plotlib import update_rc_params
import os
import numpy as np
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
os.makedirs("figures/", exist_ok=True)
import matplotlib as mpl
import matplotlib.pyplot as plt
import datetime as dt

import sys
sys.path.append("../Superstorms/py/")


def setups(size = 15):
    plt.style.use(["science", "ieee"])
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [
        "Tahoma", "DejaVu Sans",
        "Lucida Grande", "Verdana"
    ]
    mpl.rcParams.update(
        {"xtick.labelsize": size, "ytick.labelsize": size, "font.size": size}
    )
    return

setups()
bth = pd.read_csv("dataset/lat_long_bathymetry.csv")
segments = [(0,32), (32,50), (50,60), (60,170), (170,330), (330,410), (410,-1)]
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
          "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan",
          "gold", "limegreen", "darkviolet", "crimson", "teal", 
          "peru", "orchid", "slategray", "salmon", "darkkhaki"]

update_rc_params(science=True)
import matplotlib.pyplot as plt
fig, ax = plt.subplots(dpi=180, figsize=(5,3), nrows=1, ncols=1)
ax.plot(bth.index, -1*bth["bathymetry.meters"], color="k", lw=0.6)
dSegments = []
for i, seg in enumerate(segments):
    o = bth.iloc[seg[0]:seg[1]]
    ax.plot(o.index, -1*o["bathymetry.meters"], marker=".", ls="None", ms=1.2, color=colors[i])
    print(f"Initial {o.lat.round(2).iloc[0]},{o.lon.round(2).iloc[0]}")
    print(f"Final {o.lat.round(2).iloc[-1]},{o.lon.round(2).iloc[-1]}")
    dSegments.append([o.lat.round(2).iloc[0], o.lon.round(2).iloc[0]])
    #dSegments.append([o.lat.round(2).iloc[-1], o.lon.round(2).iloc[-1]])
dSegments.append([o.lat.round(2).iloc[-1], o.lon.round(2).iloc[-1]])
ax.set_xticks([])
ax.set_xlabel("")
ax.set_xlim(bth.index[0], bth.index[-1])
ax.axhline(0, ls="--", lw=0.4, color="b", alpha=0.7)
ax.set_ylim(-5000, 500)
ax.set_ylabel("Depths, m")
fig.savefig("figures/bathymetry.png")

D, H, Z = (
    pd.read_csv("dataset/D[Eskdalemuir]-rescale-HR.csv", parse_dates=["Time"]),
    pd.read_csv("dataset/H[Eskdalemuir]-rescale-HR.csv", parse_dates=["Time"]),
    pd.read_csv("dataset/Z[Eskdalemuir]-rescale-HR.csv", parse_dates=["Time"])
)
D = D[(D.Time>=dt.datetime(1958, 2, 10, 11)) & (D.Time<dt.datetime(1958, 2, 11, 8))]
H = H[(H.Time>=dt.datetime(1958, 2, 10, 11)) & (H.Time<dt.datetime(1958, 2, 11, 8))]
Z = Z[(Z.Time>=dt.datetime(1958, 2, 10, 11)) & (Z.Time<dt.datetime(1958, 2, 11, 8))]
data = pd.DataFrame()
data["Time"], data["H"], data["D"], data["Z"] = D.Time, np.copy(H.H), np.copy(D.D), np.copy(Z.Z)
data["X"], data["Y"] = (
    data.H * np.cos(np.deg2rad(data.D)),
    data.H * np.sin(np.deg2rad(data.D))
)
data = data.set_index("Time").resample("1S").asfreq().interpolate(method="cubic").reset_index()
data["F"] = np.sqrt(data.X**2+data.Y**2+data.Z**2)
data = data[data.Time<dt.datetime(1958, 2, 11, 7, 59)]
print(data.tail(60))

major_locator=mdates.HourLocator(byhour=range(0, 24, 6))
minor_locator=mdates.HourLocator(byhour=range(0, 24, 2))
formatter=DateFormatter(r"%H^{%M}")
fig, ax = plt.subplots(dpi=240, figsize=(8,3), nrows=1, ncols=1)
ax.xaxis.set_major_locator(major_locator)
ax.xaxis.set_minor_locator(minor_locator)
ax.xaxis.set_major_formatter(formatter)
ax.plot(data.Time, data.X-np.nanmedian(data.X.iloc[0:60]), color="b", lw=0.6, ls="-", label=r"$B_x$")
ax.plot(data.Time, data.Y-np.nanmedian(data.Y.iloc[0:60]), color="r", lw=0.6, ls="-", label=r"$B_y$")
ax.plot(data.Time, data.Z-np.nanmedian(data.Z.iloc[0:60]), color="g", lw=0.6, ls="-", label=r"$B_z$")
ax.legend(loc=1)
ax.set_xlim([dt.datetime(1958, 2, 10, 11), dt.datetime(1958, 2, 11, 8)])
ax.text(0.1, 1.05, "Staring at 11 UT, 10 Feb, 1958", ha="left", va="center", transform=ax.transAxes)
ax.set_ylabel("dB (Eskdalemuir), nT")
ax.set_xlabel("Time, UT")
fig.savefig("figures/EskdalemuirXY.png")

data = data.rename(columns=dict(Time="Date"))
data[["Date", "X", "Y", "Z", "F"]].to_csv("dataset/compiled.csv", float_format="%g", index=False)


from scubas.conductivity import ConductivityProfile
profiles = ConductivityProfile.compile_bined_profiles(np.array(dSegments))
for p, seg in zip(profiles, segments):
    o = bth.iloc[seg[0]:seg[1]]
    depth = o["bathymetry.meters"].max()
    p.layers[0].thickness = depth/1e3

from utils import create_from_lat_lon
from cable import SCUBASModel
from scubas.datasets import PROFILES

land50 = PROFILES.CS_E
land50.layers[0].thickness = 50
cable0, cable1, cable2 = (
    create_from_lat_lon(
        dSegments, 
        profiles, 
    ),
    create_from_lat_lon(
        dSegments, 
        profiles, 
        left_active_termination=PROFILES.LD,
        right_active_termination=PROFILES.LD
    ),
    create_from_lat_lon(
        dSegments, 
        profiles, 
        left_active_termination=PROFILES.LD,
        right_active_termination=land50
    )
)
fig_titles=[
    "No active terminations",
    "Active terminations with Land on both sides",
    "Active terminations with Land and 50-m water"
]

segment_files = [
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"], 
]

model0, model1, model2 = (
    SCUBASModel(
        cable_name="TAT-1",
        cable_structure=cable0,
        segment_files=segment_files,
    ),
    SCUBASModel(
        cable_name="TAT-1",
        cable_structure=cable1,
        segment_files=segment_files,
    ),
    SCUBASModel(
        cable_name="TAT-1",
        cable_structure=cable2,
        segment_files=segment_files,
    )
)

model0.read_stations(["ESK"], [["dataset/compiled.csv"]])
model0.initialize_TL()
model0.run_cable_segment()

model1.read_stations(["ESK"], [["dataset/compiled.csv"]])
model1.initialize_TL()
model1.run_cable_segment()

model2.read_stations(["ESK"], [["dataset/compiled.csv"]])
model2.initialize_TL()
model2.run_cable_segment()

from plots import TimeSeriesPlot
ts = TimeSeriesPlot(
    [dt.datetime(1958, 2, 10, 11), dt.datetime(1958, 2, 11, 8)],
    major_locator=major_locator,
    minor_locator=minor_locator,
    fig_title=fig_titles[0], 
    text_size=15,
    num_subplots=2,
    formatter=DateFormatter(r"%H^{%M}"),
)
ax = ts.add_voltage(model0.cable.tot_params, lw=0.7, color="k",
        ylim=[-3000, 3000], xlabel="UT since at 11 on 10 Feb, 1958")
ax.plot(model0.cable.tot_params.index, model0.cable.tot_params["V(v)"],
        color="r", ls="-", lw=0.7, alpha=0.7, zorder=1)
ts.save("figures/runs.png")
ts.close()

from plots import TimeSeriesPlot
ts = TimeSeriesPlot(
    [dt.datetime(1958, 2, 11), dt.datetime(1958, 2, 11, 4)],
    major_locator=major_locator,
    minor_locator=minor_locator,
    fig_title="", 
    text_size=15,
    num_subplots=2,
    formatter=DateFormatter(r"%H^{%M}"),
)
ax = ts.add_voltage(model0.cable.tot_params, lw=0.7, color="k",
        ylim=[-3000, 3000], xlabel="Staring at 00 UT, 11 Feb, 1958")
ax.plot(model1.cable.tot_params.index, model1.cable.tot_params["Vt(v)"],
        color="r", ls="-", lw=0.7, alpha=0.7, zorder=3, label="AT/Land",)
ax.plot(model2.cable.tot_params.index, model2.cable.tot_params["Vt(v)"],
        color="b", ls="-", lw=0.7, alpha=0.7, zorder=5, label="AT/Land,50-m Water")
ts.save("figures/compare.png")
ts.close()