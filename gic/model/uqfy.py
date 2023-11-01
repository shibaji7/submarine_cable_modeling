
"""
Uncertainity Quantification
---------------------------
Quantifiy unceratinity in the model using following informations.
We undertand that following are the primary source of unceratinity in this
model. We estimate unceratinity in the ocean water depth along the routes,
Earth Layer, and then propagate the unceratinity through the model to 
estimate.

1. Uncerayinity in Ocean depth along each section
2. Unceratinity in the Earth layer along each section
"""

import math

import numpy as np
import pyproj
from loguru import logger

from gic.model.conductivity import ConductivityProfile


class QuantifyUnceartainity(object):
    
    def __init__(
        self, 
        control, 
        input_uq = {
            "samples_per_cable_setion": 11
        }
    ):
        """
        Parametres
        ----------
        control: Control variable holding all parameres
        input_uq: Parameters for UQ
        """
        self.control = control
        self.input_uq = input_uq
        self.cf = ConductivityProfile()
        self.cable_section_profiles = {}
        self.quntify_uncertainity_in_cable_sections()
        return
    
    def quntify_uncertainity_in_cable_sections(self):
        """
        Extract uncertainities in the inputs
        """
        logger.info("Quantifying input unceratinity along cable routes")
        geodesic = pyproj.Geod(ellps="WGS84")
        samples_per_cable_setion = self.input_uq["samples_per_cable_setion"]
        R = 6378.1
        for i, cs in enumerate(self.control.cable.cable_sections):
            logger.info(f"Estimating unceratinity along cable section {cs.sec_id}")
            edge_loc = cs.edge_loc
            ini_lat, ini_lon, fin_lat, fin_lon = (
                edge_loc.ini.lat, edge_loc.ini.lon,
                edge_loc.fin.lat, edge_loc.fin.lon
            )
            fwd_azimuth, back_azimuth, distance = \
                geodesic.inv(ini_lon, ini_lat, fin_lon, fin_lat)
            distance /= 1e3
            lats, lons = [], []
            ini_lat, ini_lon, brng = (
                math.radians(ini_lat), 
                math.radians(ini_lon),
                math.radians(fwd_azimuth),
            )
            for r in np.linspace(0, distance, samples_per_cable_setion):
                lat = math.asin( math.sin(ini_lat)*math.cos(r/R) +\
                                 math.cos(ini_lat)*math.sin(r/R)*\
                                 math.cos(fwd_azimuth)
                               )
                lon = ini_lon + math.atan2(
                    math.sin(brng)*math.sin(r/R)*math.cos(ini_lat),
                    math.cos(r/R)-math.sin(ini_lat)*math.sin(lat)
                )
                lat, lon = (
                    math.degrees(lat),
                    math.degrees(lon)
                )
                lats.append(lat)
                lons.append(lon)
            binlatlons = np.array([lats, lons])
            profiles = self.cf.compile_profiles(binlatlons.T)
            self.cable_section_profiles[cs.sec_id] = dict(
                lats=lats,
                lons=lons,
                profiles=profiles
            )
        return
    
    def cable_runs(self):
        """
        Running the cable sections based on various input stats
        """
        samples_per_cable_setion = self.input_uq["samples_per_cable_setion"]
        for sample in range(samples_per_cable_setion):
            
        return
    