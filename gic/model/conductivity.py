"""conductiviy.py: Module is used to implement Earth conductivity methods"""

__author__ = "Murphy, B.; Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import json

# import sys
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import siunits as u
from loguru import logger
from pyproj import Geod
from scipy.interpolate import RegularGridInterpolator
from scipy.io import netcdf_file

# sys.path.extend(["py/", "config/"])


def get_value(x):
    """
    Extract value from <v,u> data.
    """
    return x.value


def get_static_profile_names(fname="config/static_conductivity_profiles.json"):
    """
    Load and return profile names
    """
    with open(fname, "r") as f:
        o = json.load(f)
    profiles = []
    for p in o["profiles"].keys():
        profiles.append(
            {
                "code": p,
                "name": o["profiles"][p]["name"],
                "desciption": o["profiles"][p]["desciption"],
            }
        )
    return profiles


def fetch_static_conductivity_profiles(
    profile="OM", fname="config/static_conductivity_profiles.json"
):
    """
    Fetch static conductivity profile
    """
    with open(fname, "r") as f:
        o = json.load(f)
    t, c = o["profiles"][profile]["thickness"], o["profiles"][profile]["conductivity"]
    rf = {}
    rf["thickness"], rf["resistivity"] = np.array(t), 1.0 / np.array(c)
    return rf


def load_conductivity_params(fname="config/conductivity.json"):
    """
    Load conductivity properties
    """
    with open(fname, "r") as f:
        o = json.load(f, object_hook=lambda d: SimpleNamespace(**d))
    return o


class ConductivityProfile(object):
    """
    Class is dedicated to create conductivity profiles
    from LITHO1.0 model.

    1.  LITHO1.0 as the basis for constructing the resistivity profiles
        http://ds.iris.edu/ds/products/emc-litho10/
    """

    def __init__(self):
        self.cprop = load_conductivity_params()
        self.earth_model = self.cprop.earth_model

        unit = self.cprop.units.resistivity
        self.ohm_m = (
            u.__dict__[self.cprop.units.resistivity.split("-")[0]]
            * u.__dict__[self.cprop.units.resistivity.split("-")[1]]
        )
        self.km = u.DerivedUnit("kilometer", "km", [(u.m, 1000)], "length")

        # NOTE these are specified in terms of resistivity, in ohm-m
        self.seawater_resistivity = self.cprop.seawater_resistivity
        self.sediment_resistivity = self.cprop.sediment_resistivity
        self.crust_resistivity = self.cprop.crust_resistivity
        self.lithosphere_resistivity = self.cprop.lithosphere_resistivity
        self.asthenosphere_resistivity = self.cprop.asthenosphere_resistivity
        self.transition_zone_resistivity = self.cprop.transition_zone_resistivity
        self.lower_mantle_resistivity = self.cprop.lower_mantle_resistivity

        # fixed layer depths (below the lithosphere), in km
        self.transition_zone_top = self.cprop.transition_zone_top
        self.transition_zone_bot = self.cprop.transition_zone_bot
        self.profile_max_depth = self.cprop.profile_max_depth

        # things that control how this code works in detail...
        # Options: "nearest" or "linear"
        self.grid_interpolation_method = self.cprop.grid_interpolation_method

        # Load netcdf file
        self.load_earth_model()
        return

    def load_earth_model(self):
        """
        Returns a dict containing the pertinent parts of the Earth model
        Dict entries are callable interpolation objects
        """
        filename = "config/" + self.earth_model
        if not os.path.exists(filename):
            uri = self.cprop.uri
            os.system(f"wget -O {filename} {uri}")
        with netcdf_file(filename) as f:
            latitude = np.copy(f.variables["latitude"][:])
            longitude = np.copy(f.variables["longitude"][:])
            # base of lithosphere (top of asthenosphere)
            asthenospheric_mantle_top_depth = np.copy(
                f.variables["asthenospheric_mantle_top_depth"][:]
            )
            # top/bottom of mantle lithosphere
            lithosphere_bottom_depth = np.copy(f.variables["lid_bottom_depth"][:])
            lithosphere_top_depth = np.copy(f.variables["lid_top_depth"][:])
            # crustal layers
            lower_crust_bottom_depth = np.copy(
                f.variables["lower_crust_bottom_depth"][:]
            )
            lower_crust_top_depth = np.copy(f.variables["lower_crust_top_depth"][:])
            middle_crust_bottom_depth = np.copy(
                f.variables["middle_crust_bottom_depth"][:]
            )
            middle_crust_top_depth = np.copy(f.variables["middle_crust_top_depth"][:])
            upper_crust_bottom_depth = np.copy(
                f.variables["upper_crust_bottom_depth"][:]
            )
            upper_crust_top_depth = np.copy(f.variables["upper_crust_top_depth"][:])
            # sediment layers
            lower_sediments_bottom_depth = np.copy(
                f.variables["lower_sediments_bottom_depth"][:]
            )
            lower_sediments_top_depth = np.copy(
                f.variables["lower_sediments_top_depth"][:]
            )
            middle_sediments_bottom_depth = np.copy(
                f.variables["middle_sediments_bottom_depth"][:]
            )
            middle_sediments_top_depth = np.copy(
                f.variables["middle_sediments_top_depth"][:]
            )
            upper_sediments_bottom_depth = np.copy(
                f.variables["upper_sediments_bottom_depth"][:]
            )
            upper_sediments_top_depth = np.copy(
                f.variables["upper_sediments_top_depth"][:]
            )
            # water levels
            water_bottom_depth = np.copy(f.variables["water_bottom_depth"][:])
            water_top_depth = np.copy(f.variables["water_top_depth"][:])

        self.lithosphere_model = {
            "latitude": latitude,
            "longitude": longitude,
            "asthenospheric_mantle_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                asthenospheric_mantle_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lithosphere_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                lithosphere_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lithosphere_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                lithosphere_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lower_crust_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                lower_crust_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lower_crust_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                lower_crust_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "middle_crust_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                middle_crust_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "middle_crust_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                middle_crust_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "upper_crust_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                upper_crust_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "upper_crust_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                upper_crust_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lower_sediments_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                lower_sediments_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lower_sediments_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                lower_sediments_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "middle_sediments_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                middle_sediments_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "middle_sediments_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                middle_sediments_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "upper_sediments_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                upper_sediments_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "upper_sediments_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                upper_sediments_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "water_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                water_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "water_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                water_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
        }
        return

    def get_interpolation_points(self, pt0, pt1):
        globe = Geod(ellps="WGS84")

        # somewhat janky way of figuring out how many interp points to have between
        # the bin edges... basic idea is that we want one point spaced every ~1 deg,
        # since that"s the resolution of LITHO1.0
        npts = int(round(np.sqrt((pt1[1] - pt0[1]) ** 2.0 + (pt1[0] - pt0[0]) ** 2.0)))
        if npts < 1:
            npts = 1

        # obtain equally spaced points along a geodesic (doesn"t include end points)
        latlons = globe.npts(pt0[1], pt0[0], pt1[1], pt1[0], npts)

        # add the bin edges to the lat/lon list... note this array is now [lon, lat]
        interpolation_points = np.vstack(
            (np.array([pt0[1], pt0[0]]), np.array(latlons), np.array([pt1[1], pt1[0]]))
        )

        # now swap back to [lat, lon]
        interpolation_points = np.vstack(
            (interpolation_points[:, 1], interpolation_points[:, 0])
        ).T

        return interpolation_points

    def get_water_layer(self, lithosphere_model, pts):
        """
        Interpolates and averages seawater layer thickness within the given bin
        """

        water_top_values = lithosphere_model["water_top_depth"](pts)
        water_bot_values = lithosphere_model["water_bottom_depth"](pts)

        water_top = np.nanmean(water_top_values)
        water_bot = np.nanmean(water_bot_values)
        water_thk = water_bot - water_top

        # sanity check... water_top can be slightly different than zero,
        # but shouldn"t be that much different
        if np.absolute(water_top) > 0.01:
            logger.warning(f"PROBLEM: water_top doesn't make sense: {water_top}")

        return water_thk

    def get_sediment_layer(self, lithosphere_model, pts):
        """
        Interpolates and averages sediment layer thickness within the given bin
        """
        upper_sed_top_values = lithosphere_model["upper_sediments_top_depth"](pts)
        upper_sed_bot_values = lithosphere_model["upper_sediments_bottom_depth"](pts)
        upper_sed_thk_values = upper_sed_bot_values - upper_sed_top_values

        middle_sed_top_values = lithosphere_model["middle_sediments_top_depth"](pts)
        middle_sed_bot_values = lithosphere_model["middle_sediments_bottom_depth"](pts)
        middle_sed_thk_values = middle_sed_bot_values - middle_sed_top_values

        lower_sed_top_values = lithosphere_model["lower_sediments_top_depth"](pts)
        lower_sed_bot_values = lithosphere_model["lower_sediments_bottom_depth"](pts)
        lower_sed_thk_values = lower_sed_bot_values - lower_sed_top_values

        sed_thk_values = np.nansum(
            np.vstack(
                (upper_sed_thk_values, middle_sed_thk_values, lower_sed_thk_values)
            ),
            axis=0,
        )

        sed_thk = np.nanmean(sed_thk_values)

        return sed_thk

    def get_crust_layer(self, lithosphere_model, pts):
        """
        Interpolates and averages crust layer thickness within the given bin
        """
        upper_crust_top_values = lithosphere_model["upper_crust_top_depth"](pts)
        upper_crust_bot_values = lithosphere_model["upper_crust_bottom_depth"](pts)
        upper_crust_thk_values = upper_crust_bot_values - upper_crust_top_values

        middle_crust_top_values = lithosphere_model["middle_crust_top_depth"](pts)
        middle_crust_bot_values = lithosphere_model["middle_crust_bottom_depth"](pts)
        middle_crust_thk_values = middle_crust_bot_values - middle_crust_top_values

        lower_crust_top_values = lithosphere_model["lower_crust_top_depth"](pts)
        lower_crust_bot_values = lithosphere_model["lower_crust_bottom_depth"](pts)
        lower_crust_thk_values = lower_crust_bot_values - lower_crust_top_values

        crust_thk_values = np.nansum(
            np.vstack(
                (
                    upper_crust_thk_values,
                    middle_crust_thk_values,
                    lower_crust_thk_values,
                )
            ),
            axis=0,
        )

        crust_thk = np.nanmean(crust_thk_values)

        return crust_thk

    def get_lithosphere_layer(self, lithosphere_model, pts):
        """
        Interpolates and averages mantle lithosphere thickness within the given bin
        """

        lithosphere_top_values = lithosphere_model["lithosphere_top_depth"](pts)
        lithosphere_bot_values = lithosphere_model["lithosphere_bottom_depth"](pts)
        asthenosphere_top_values = lithosphere_model["asthenospheric_mantle_top_depth"](
            pts
        )

        litho_top = np.nanmean(lithosphere_top_values)
        litho_bot = np.nanmean(lithosphere_bot_values)
        litho_thk = litho_bot - litho_top
        astheno_top = np.nanmean(asthenosphere_top_values)

        # sanity check... make sure top of asthenosphere is the same as bottom of the lithosphere
        if litho_bot != astheno_top:
            logger.warning(
                f"PROBLEM: lithosphere-asthenosphere dont line up: {litho_bot}, {astheno_top}"
            )

        return litho_thk

    def get_upper_mantle_layer(self, lithosphere_model, pts):
        """
        Interpolates and averages upper mantle thickness within the given bin
        """
        asthenosphere_top_values = lithosphere_model["asthenospheric_mantle_top_depth"](
            pts
        )
        astheno_top = np.nanmean(asthenosphere_top_values)
        astheno_thk = self.transition_zone_top - astheno_top
        return astheno_thk

    def get_transition_zone_layer(self):
        """
        Return fixed mantle transition zone thickness
        """
        return self.transition_zone_bot - self.transition_zone_top

    def get_lower_mantle_layer(self):
        """
        Return fixed lower mantle thickness
        """
        return self.profile_max_depth - self.transition_zone_bot

    def _compile_profile_(self, pts):
        """
        Compile profile for a list of points
        """
        # now progress through the model layers from near-surface to deep Earth...
        # first the water layer
        water_thk = self.get_water_layer(self.lithosphere_model, pts)
        # sediment layer
        sed_thk = self.get_sediment_layer(self.lithosphere_model, pts)
        # crust layer
        crust_thk = self.get_crust_layer(self.lithosphere_model, pts)
        # mantle lithosphere layer
        litho_thk = self.get_lithosphere_layer(self.lithosphere_model, pts)
        # asthenosphere layer
        astheno_thk = self.get_upper_mantle_layer(self.lithosphere_model, pts)
        # transition zone layer
        tz_thk = self.get_transition_zone_layer()
        # lower mantle layer
        lm_thk = self.get_lower_mantle_layer()
        resistivity_profile = np.array(
            [
                [water_thk, self.seawater_resistivity],
                [sed_thk, self.sediment_resistivity],
                [crust_thk, self.crust_resistivity],
                [litho_thk, self.lithosphere_resistivity],
                [astheno_thk, self.asthenosphere_resistivity],
                [tz_thk, self.transition_zone_resistivity],
                [lm_thk, self.lower_mantle_resistivity],
            ]
        )
        rf = pd.DataFrame()
        rf["thickness"], rf["resistivity"] = (
            resistivity_profile[:, 0] * self.km,
            resistivity_profile[:, 1] * self.ohm_m,
        )
        return rf

    def compile_profiles(self, latlons, kind="rounded"):
        """
        Compile profiles for a set of latlons
        """
        profiles = []
        for latlon in latlons:
            if kind == "rounded":
                latlon = np.rint(latlon)
            logger.info(f"Lat-lon: {latlon}")
            rf = self._compile_profile_(latlon)
            logger.info(f"Compiled Profile \n {rf}")
            profiles.append(rf)
        return profiles

    def compile_bin_profiles(self, binlatlons):
        """
        Compile profiles for a set of binned latlons
        """
        profiles = []
        nbins = len(binlatlons) - 1
        for i in range(nbins):
            ipts = self.get_interpolation_points(binlatlons[i, :], binlatlons[i + 1, :])
            rf = self._compile_profile_(ipts)
            logger.info(f"Compiled Profile \n {rf}")
            profiles.append(rf)
        return profiles


if __name__ == "__main__":
    # Testing codes
    latlons = np.array([[48, -167], [-36, 160], [54, 3]])
    cf = ConductivityProfile()
    cf.compile_profiles(latlons)
