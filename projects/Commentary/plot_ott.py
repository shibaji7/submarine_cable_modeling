import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import matplotlib.dates as mdates
import datetime as dt


import numpy as np
import pandas as pd
import glob

files = glob.glob("datasets/*.csv")
files.sort()
datasets = pd.concat([pd.read_csv(f, parse_dates=["datetime"]) for f in files])
datasets.x = datasets.x-np.mean(datasets.x.iloc[:60])
datasets.y = datasets.y-np.mean(datasets.y.iloc[:60])
datasets.z = datasets.z-np.mean(datasets.z.iloc[:60])

fig = plt.figure(dpi=300, figsize=(6,1.5))
ax = fig.add_subplot(111)
date_format = mdates.DateFormatter("%H")
ax.xaxis.set_major_formatter(date_format)
ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
ax.set_xlim(dt.datetime(1989,3,12), dt.datetime(1989,3,15))
ax.set_ylim(-2000, 2000)
ax.set_ylabel(r"B, nT / OTT")
ax.set_xlabel(r"Time, UT since 12 March 1989")
ax.plot(datasets.datetime, datasets.x, ls="-", lw=0.6, color="r", label=r"$B_x$")
ax.plot(datasets.datetime, datasets.y, ls="-", lw=0.6, color="k", label=r"$B_y$")
ax.plot(datasets.datetime, datasets.z, ls="-", lw=0.6, color="b", label=r"$B_z$")
ax.legend(loc=1)
fig.savefig("datasets.png", bbox_inches="tight")