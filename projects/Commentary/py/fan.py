#!/usr/bin/env python

"""
    fanUtils.py: module to plot Fan plots with various transformation
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]

import cartopy
import matplotlib.ticker as mticker
from cartobase import CartoBase, setsize
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER


class CartoDataOverlay(object):

    def __init__(
        self,
        date,
        fig_title=None,
        nrows=3,
        ncols=3,
        coord="geo",
        cb=True,
        central_longitude=-30.0,
        central_latitude=45.0,
        extent=[-130, -60, 20, 80],
        plt_lats=np.arange(20, 80, 10),
        txt_size=8,
    ):
        setsize(txt_size)
        self.cb = cb
        self.date = date
        self.nrows, self.ncols = nrows, ncols
        self._num_subplots_created = 0
        self.fig = plt.figure(figsize=(3 * ncols, 3 * nrows), dpi=300)
        self.coord = coord
        self.central_longitude = central_longitude
        self.central_latitude = central_latitude
        self.plt_lats = plt_lats
        self.extent = extent
        self.fig_title = fig_title
        return

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight", facecolor=(1, 1, 1, 1))
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return

    def date_string(self, label_style="web"):
        # Set the date and time formats
        dfmt = "%d %b %Y" if label_style == "web" else "%d %b %Y,"
        tfmt = "%H:%M"
        stime = self.date
        date_str = "{:{dd} {tt}} UT".format(stime, dd=dfmt, tt=tfmt)
        return date_str

    def add_axes(self):
        """
        Instatitate figure and axes labels
        """
        self._num_subplots_created += 1
        proj = cartopy.crs.PlateCarree(
            central_longitude=self.central_longitude,
        )
        ax = self.fig.add_subplot(
            100 * self.nrows + 10 * self.ncols + self._num_subplots_created,
            projection="CartoBase",
            map_projection=proj,
            coords=self.coord,
            plot_date=self.date,
        )
        ax.overaly_coast_lakes(lw=0.4, alpha=0.4)
        plt_lons = np.arange(-180, 181, 5)
        mark_lons = np.arange(self.extent[0], self.extent[1], 5)
        plt_lats = self.plt_lats
        ax.set_extent(self.extent, crs=cartopy.crs.PlateCarree())
        gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), linewidth=0.2)
        gl.xlocator = mticker.FixedLocator(plt_lons)
        gl.ylocator = mticker.FixedLocator(plt_lats)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.n_steps = 90
        ax.mark_latitudes(plt_lats, fontsize="xx-small", color="k", lon_location=-90)
        ax.mark_longitudes(mark_lons, fontsize="xx-small", color="k")
        self.proj = proj
        self.geo = cartopy.crs.PlateCarree()
        if self._num_subplots_created == 1:
            ax.text(
                -0.02,
                0.99,
                "Coord: Geo",
                ha="center",
                va="top",
                transform=ax.transAxes,
                fontsize="x-small",
                rotation=90,
            )
            ax.text(
                0.05,
                1.05,
                (
                    f"{self.date_string()} / {self.fig_title}"
                    if self.fig_title
                    else f"{self.date_string()}"
                ),
                ha="left",
                va="center",
                transform=ax.transAxes,
                fontsize="x-small",
            )
        return ax
    
def get_bthmetry():
    # from scubas.coductivity import ConductivityProfile
    # cp = ConductivityProfile()
    from scipy.io import netcdf_file
    filename = ".scubas_config/LITHO1.0.nc"
    with netcdf_file(filename) as f:
        latitude = np.copy(f.variables["latitude"][:])
        longitude = np.copy(f.variables["longitude"][:])

        # water levels
        water_bottom_depth = np.copy(f.variables["water_bottom_depth"][:])
        water_top_depth = np.copy(f.variables["water_top_depth"][:])
        dwater = water_bottom_depth-water_top_depth
        dwater[dwater<=0] = np.nan
        # dwater = np.ma.masked_invalid(dwater)
    return (latitude, longitude, dwater)

if __name__ == "__main__":
    cb = CartoDataOverlay(date=dt.datetime(1958,2,11))
    ax = cb.add_axes()
    bth = pd.read_csv("dataset/lat_long_bathymetry.csv")
    xyz = cb.proj.transform_points(
        cb.geo, bth.lon, bth.lat
    )
    ax.plot(
        xyz[:, 0], xyz[:, 1],
        ls="-", lw=0.8, color="r",
        transform=cb.proj
    )
    (latitude, longitude, dwater) = get_bthmetry()
    lats, lons = np.meshgrid(latitude, longitude)
    xyz = cb.proj.transform_points(
        cb.geo, lons, lats
    )
    im = ax.pcolormesh(
        xyz[:,:,0], xyz[:,:,1],
        dwater.T, cmap="Blues",
        transform=cb.proj, vmax=4,
        vmin=0
    )    
    cpos = [1.04, 0.1, 0.025, 0.8]
    cax = ax.inset_axes(cpos, transform=ax.transAxes)
    cbr = cb.fig.colorbar(im, ax=ax, cax=cax)
    cbr.set_label("Water Depth (km)", color="k", fontsize="x-small")
    cbr.set_ticks(np.linspace(0, 4, 5))
    cbr.ax.tick_params(which="both", colors="k")
    cbr.outline.set_edgecolor("k")
    cb.save("figures/routes.png")
    cb.close()
    pass