import pandas as pd
import numpy as np
import datetime as dt
import glob
from scipy import constants as C

from scubas.datasets import PROFILES
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
        ):
        self.cable_name = cable_name
        self.cable_structure = cable_structure
        self.initialize_TL()
        self.run_cable_segment()
        return
    
    def initialize_TL(self):
        self.tlines = []
        for seg in self.cable_structure.cable_seg:
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
                ).compile_oml(FRD_files),
            )
        return
    
    def run_cable_segment(self):
        # Running the cable operation
        self.cable = Cable(self.tlines, self.tlines[0].components)
        return
    
    def plot_TS_with_others(self):
        return