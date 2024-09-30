import pandas as pd
import datetime as dt
import numpy as np

from scubas.datasets import PROFILES

from event_runs import (
    read_bfield_frames, 
    compile_cable_to_calculate_parameters,
    initialize_mag_fold
)
from event_plots import plot_potential


initialize_mag_fold()
FRD, HAD, STJ = (
    [".scubas_config/frd20240510psec.sec.txt"],
    [".scubas_config/had20240510psec.sec.txt"], 
    [".scubas_config/stj20240510vsec.sec.txt"]
)
stns, dates = (
    ["FRD", "STJ", "HAD"], 
    [dt.datetime(2024,5,10,12), dt.datetime(2024,5,12)]
)
frames = read_bfield_frames(stns, [FRD, HAD, STJ])

profiles = [
    PROFILES.CS_W, PROFILES.DO_1, PROFILES.DO_2, PROFILES.DO_3,
    PROFILES.DO_4, PROFILES.DO_5, PROFILES.MAR, PROFILES.DO_6,
    PROFILES.CS_E
]
tlines, cable = compile_cable_to_calculate_parameters(FRD, STJ, HAD, profiles)
plot_potential(
    frames, dates,
    cable.tot_params.index, [cable.tot_params["V(v)"]]
)
# Run for initial spike event 17:10 UT 10 May
# Run for Sharp turning event 2:30 UT 11 May
# Run for sharp turning event 12 UT 11 May

# Plot of collrelation analysis among stations
# Plot DMSP data-plots and SuperMag dataset for infernece
# Email MoM to group