import sys
sys.path.append("py/*")
import utils

import numpy as np
from scubas.utils import fft, ifft
import matplotlib as mpl
import matplotlib.dates as mdates
import datetime as dt
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(["science", "ieee"])
plt.rcParams.update({
    "text.usetex": False,
})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
plt.rcParams.update({
    "font.size": 10
})

ds = utils.get_benchmark_datasets()
ds_b = utils.get_nerc_dataset()
select_rows = 1
ds = ds.iloc[::select_rows]
site = utils.sites[4]
print(ds.head(), len(ds))
print(ds_b.head(), len(ds_b))

Bxf, f = fft(ds.x, 10*select_rows)
Ets = dict(
    Y=-ifft(np.array(utils.calcTFx(site, freqs=f).E2B) * Bxf),
)
Byf, f = fft(ds.y, 10*select_rows)
Ets["X"] = ifft(np.array(utils.calcTFx(site, freqs=f).E2B) * Byf)

fig, start = plt.figure(dpi=300, figsize=(5,7)), 311
ax = fig.add_subplot(start)
date_format = mdates.DateFormatter("%H")
ax.xaxis.set_major_formatter(date_format)
ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,3))
ax.set_ylim(-4000, 2000)
ax.set_ylabel(r"B, nT")
ax.plot(ds.datetime, ds.x, ls="-", lw=0.6, color="r", label="$B_x$")
ax.plot(ds.datetime, ds.y, ls="-", lw=0.6, color="b", label="$B_y$")
ax.legend(loc=1)
ax.set_xticklabels([])
ax.text(0.05, 0.95, "(a)", ha="left", va="top", transform=ax.transAxes)

ax = fig.add_subplot(start+1)
date_format = mdates.DateFormatter("%H")
ax.xaxis.set_major_formatter(date_format)
ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
ax.set_xlim(dt.datetime(1989,3,13), dt.datetime(1989,3,14,3))
ax.set_ylim(-8, 8)
ax.set_ylabel(r"E, V/km")
ax.plot(ds.datetime, Ets["X"]/1e3, ls="-", lw=0.6, color="r")
ax.plot(ds.datetime, Ets["Y"]/1e3, ls="-", lw=0.6, color="b")
# ax.plot(ds.datetime, utils.dBdt_fft(ds.y, sqrt=True), ls="-", lw=0.6, color="k")
ax.text(0.05, 0.95, "(b) Our simulation", ha="left", va="top", transform=ax.transAxes)
ax.set_xlabel(r"Time, UT since 0 UT on 13 March 1989")


ax = fig.add_subplot(start+2)
ax.plot(ds_b.Time/(3600), ds_b["E_NS(mV/km)"]/1e3, ls="-", lw=0.6, color="r")
ax.plot(ds_b.Time/(3600), ds_b["E_EW(mV/km)"]/1e3, ls="-", lw=0.6, color="b")
ax.set_xlim(0, 27)
ax.set_ylim(-8, 8)
ax.set_ylabel(r"E, V/km")
ax.set_xlabel(r"Time, Hours since 0 UT on 13 March 1989")
ax.text(0.05, 0.95, "(c) NERC simulated", ha="left", va="top", transform=ax.transAxes)

fig.subplots_adjust(wspace=0.6, hspace=0.5)
fig.savefig("figures/Figure04a.png", bbox_inches="tight")
# print(np.abs(Ets["X"]/1e3).max())
