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

<<<<<<< HEAD
=======
import pandas as pd

>>>>>>> a8429f5bf3936f33668b28c9a7eedbb37d96a9e7

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
<<<<<<< HEAD
        extent=[-90, 10, 40, 90],
        plt_lats=np.arange(40, 90, 10),
=======
        extent=[-90, 10, 20, 80],
        plt_lats=np.arange(20, 80, 10),
>>>>>>> a8429f5bf3936f33668b28c9a7eedbb37d96a9e7
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
<<<<<<< HEAD
        proj = cartopy.crs.NorthPolarStereo(central_longitude=self.central_longitude)
=======
        proj = cartopy.crs.PlateCarree(
            central_longitude=self.central_longitude,
        )
>>>>>>> a8429f5bf3936f33668b28c9a7eedbb37d96a9e7
        ax = self.fig.add_subplot(
            100 * self.nrows + 10 * self.ncols + self._num_subplots_created,
            projection="CartoBase",
            map_projection=proj,
            coords=self.coord,
            plot_date=self.date,
        )
        ax.overaly_coast_lakes(lw=0.4, alpha=0.4)
<<<<<<< HEAD
        plt_lons = np.arange(-180, 181, 15)
        mark_lons = np.arange(self.extent[0], self.extent[1], 15)
=======
        plt_lons = np.arange(-180, 181, 20)
        mark_lons = np.arange(self.extent[0], self.extent[1], 20)
>>>>>>> a8429f5bf3936f33668b28c9a7eedbb37d96a9e7
        plt_lats = self.plt_lats
        ax.set_extent(self.extent, crs=cartopy.crs.PlateCarree())
        gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), linewidth=0.2)
        gl.xlocator = mticker.FixedLocator(plt_lons)
        gl.ylocator = mticker.FixedLocator(plt_lats)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.n_steps = 90
        ax.mark_latitudes(plt_lats, fontsize="xx-small", color="k")
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
    
<<<<<<< HEAD

from glob import glob

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader


def load_bathymetry(zip_file_url):
    """Read zip file from Natural Earth containing bathymetry shapefiles"""
    # Download and extract shapefiles
    import io
    import zipfile

    import requests
    r = requests.get(zip_file_url)
    z = zipfile.ZipFile(io.BytesIO(r.content))
    z.extractall("dataset/ne_10m_bathymetry_all/")

    # Read shapefiles, sorted by depth
    shp_dict = {}
    files = glob('dataset/ne_10m_bathymetry_all/*.shp')
    assert len(files) > 0
    files.sort()
    depths = []
    for f in files:
        depth = '-' + f.split('_')[-1].split('.')[0]  # depth from file name
        depths.append(depth)
        bbox = (90, -15, 160, 60)  # (x0, y0, x1, y1)
        nei = shpreader.Reader(f, bbox=bbox)
        shp_dict[depth] = nei
    depths = np.array(depths)[::-1]  # sort from surface to bottom
    return depths, shp_dict
    
if __name__ == "__main__":
    depths_str, shp_dict = load_bathymetry(
        'https://naturalearth.s3.amazonaws.com/' +
        '10m_physical/ne_10m_bathymetry_all.zip')
    depths = depths_str.astype(int)

    N = len(depths)
    nudge = 0.01  # shift bin edge slightly to include data
    boundaries = [min(depths)] + sorted(depths+nudge)  # low to high
    norm = matplotlib.colors.BoundaryNorm(boundaries, N)
    blues_cm = matplotlib.colormaps['Blues_r'].resampled(N)
    colors_depths = blues_cm(norm(depths))


    import pandas as pd
    data = pd.read_csv("dataset/lat_long_bathymetry.csv")
    
    # Construct a discrete colormap with colors corresponding to each depth
    depths = depths_str.astype(int)
    N = len(depths)
    nudge = 0.01  # shift bin edge slightly to include data
    boundaries = [min(depths)] + sorted(depths+nudge)  # low to high
    norm = matplotlib.colors.BoundaryNorm(boundaries, N)
    blues_cm = matplotlib.colormaps['Blues_r'].resampled(N)
    colors_depths = blues_cm(norm(depths))

    # Set up plot
    subplot_kw = {'projection': ccrs.LambertCylindrical()}
    fig, ax = plt.subplots(subplot_kw=subplot_kw, figsize=(9, 7))
    ax.set_extent([-90, 160, -15, 90], crs=ccrs.PlateCarree())  # x0, x1, y0, y1

    # Iterate and plot feature for each depth level
    for i, depth_str in enumerate(depths_str):
        ax.add_geometries(shp_dict[depth_str].geometries(),
                          crs=ccrs.PlateCarree(),
                          color=colors_depths[i])

    # Add standard features
    ax.add_feature(cfeature.LAND, color='grey')
    ax.coastlines(lw=1, resolution='110m')
    ax.gridlines(draw_labels=False)
    ax.set_position([0.03, 0.05, 0.8, 0.9])

    # Add custom colorbar
    axi = fig.add_axes([0.85, 0.1, 0.025, 0.8])
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    sm = plt.cm.ScalarMappable(cmap=blues_cm, norm=norm)
    fig.colorbar(mappable=sm,
                 cax=axi,
                 spacing='proportional',
                 extend='min',
                 ticks=depths,
                 label='Depth (m)')

    # Convert vector bathymetries to raster (saves a lot of disk space)
    # while leaving labels as vectors
    ax.set_rasterized(True)
    fig.savefig("dataset/routes.png")
=======
def get_bthmetry():
    # from scubas.conductivity import ConductivityProfile
    # cp = ConductivityProfile()
    # latitude, longitude = (
    #     np.linspace(-90, 90, 91),
    #     np.linspace(-180, 180, 181)
    # )
    # pts = np.array([[latitude[0], longitude[0]]])
    # prof = cp.compile_profile(pts)
    # print(prof)
    # water_bot_values = cp.lithosphere_model["water_bottom_depth"](pts)
    from scipy.io import netcdf_file
    filename = ".scubas_config/LITHO1.0.nc"
    from scipy.interpolate import RegularGridInterpolator
    with netcdf_file(filename) as f:
        latitude = np.copy(f.variables["latitude"][:])
        longitude = np.copy(f.variables["longitude"][:])

        # water levels
        water_bottom_depth = np.copy(f.variables["water_bottom_depth"][:])
        water_top_depth = np.copy(f.variables["water_top_depth"][:])
        dwater = water_bottom_depth-water_top_depth
        dwater[dwater<=0] = np.nan
        dwater = np.ma.masked_invalid(dwater)
        xg, yg = np.meshgrid(longitude, latitude, indexing='ij', sparse=True)
        interp = RegularGridInterpolator((xg, yg), dwater)
    return (latitude, longitude, interp(np.array([longitude, latitude])))

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
    (latitude, longitude, water) = get_bthmetry()
    lats, lons = np.meshgrid(latitude, longitude)
    # xyz = cb.proj.transform_points(
    #     cb.geo, lons, lats
    # )
    # im = ax.pcolormesh(
    #     xyz[:,:,0], xyz[:,:,1],
    #     dwater.T, cmap="Blues",
    #     transform=cb.proj, vmax=4,
    #     vmin=0
    # )
    # cpos = [1.04, 0.1, 0.025, 0.8]
    # cax = ax.inset_axes(cpos, transform=ax.transAxes)
    # cbr = cb.fig.colorbar(im, ax=ax, cax=cax)
    # cbr.set_label("Water Depth (km)", color="k", fontsize="x-small")
    # cbr.set_ticks(np.linspace(0, 4, 5))
    # cbr.ax.tick_params(which="both", colors="k")
    # cbr.outline.set_edgecolor("k")
    cb.save("figures/routes.png")
    cb.close()
>>>>>>> a8429f5bf3936f33668b28c9a7eedbb37d96a9e7
