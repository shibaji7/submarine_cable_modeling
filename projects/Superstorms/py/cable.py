import pandas as pd
import numpy as np
import datetime as dt
import glob
from scipy import constants as C
from loguru import logger

from scubas.cables import TransmissionLine, Cable
import os
import sys
sys.path.append("py/")

from utils import get_cable_informations

class SCUBASModel(object):
    
    def __init__(
            self, 
            cable_name="TAT-8",
            cable_structure=get_cable_informations(),
            segment_files=[],
        ):
        self.cable_name = cable_name
        self.cable_structure = cable_structure
        self.segment_files = segment_files
        logger.info(f"Initialize {cable_name}")
        return
    
    def initialize_TL(self):
        self.tlines = []
        for i, seg in enumerate(self.cable_structure.cable_seg):
            self.tlines.append(
                TransmissionLine(
                    sec_id=seg["sec_id"],
                    directed_length=dict(
                        edge_locations=dict(
                            initial=seg["initial"], 
                            final=seg["final"]
                        )
                    ),
                    elec_params=dict(
                        site=seg["site"],
                        width=seg["width"],
                        flim=seg["flim"],
                    ),
                    active_termination=seg["active_termination"],
                ).compile_oml(self.segment_files[i]),
            )
            print(self.tlines[i].model)
            print(self.tlines[i].model.Efield.head())
        return
    
    def run_cable_segment(self):
        # Running the cable operation
        logger.info(f"Components: {self.tlines[0].components}")
        self.cable = Cable(self.tlines, self.tlines[0].components)
        return
    
    def plot_TS_with_others(self):
        return