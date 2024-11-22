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
D = D[(D.Time>=dt.datetime(1958, 2, 10, 11)) & (D.Time<=dt.datetime(1958, 2, 11, 8))]
H = H[(H.Time>=dt.datetime(1958, 2, 10, 11)) & (H.Time<=dt.datetime(1958, 2, 11, 8))]
Z = Z[(Z.Time>=dt.datetime(1958, 2, 10, 11)) & (Z.Time<=dt.datetime(1958, 2, 11, 8))]
data = pd.DataFrame()
data["Time"], data["H"], data["D"], data["Z"] = D.Time, np.copy(H.H), np.copy(D.D), np.copy(Z.Z)
data["X"], data["Y"] = (
    data.H * np.cos(np.deg2rad(data.D)),
    data.H * np.sin(np.deg2rad(data.D))
)
print(data.head(70))
data = data.set_index("Time").resample("1S").asfreq().interpolate(method="cubic").reset_index()
data["F"] = np.sqrt(data.X**2+data.Y**2+data.Z**2)

major_locator=mdates.HourLocator(byhour=range(0, 24, 12))
minor_locator=mdates.HourLocator(byhour=range(0, 24, 1))
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
cable = create_from_lat_lon(dSegments, profiles)

segment_files = [
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"],
    ["dataset/compiled.csv"], 
]

model = SCUBASModel(
    cable_name="TAT-1",
    cable_structure=cable,
    segment_files=segment_files,
)
model.read_stations(["ESK"], [["dataset/compiled.csv"]])
model.initialize_TL()