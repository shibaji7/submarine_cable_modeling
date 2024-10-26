import pandas as pd
import datetime as dt
import numpy as np
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

from scubas.datasets import PROFILES

from event_runs import (
    read_bfield_frames, 
    compile_cable_to_calculate_parameters,
    initialize_mag_fold,
    _load_omni_
)
from event_plots import TimeSeriesPlot


initialize_mag_fold()
FRD, HAD, STJ = (
    [
        ".scubas_config/1989/FRD_19890312_XYZ.txt",
        ".scubas_config/1989/FRD_19890313_XYZ.txt",
        ".scubas_config/1989/FRD_19890314_XYZ.txt"
    ],
    [
        ".scubas_config/1989/HAD_19890312_HDZ.txt",
        ".scubas_config/1989/HAD_19890313_HDZ.txt",
        ".scubas_config/1989/HAD_19890314_HDZ.txt" 
    ], 
    [
        ".scubas_config/1989/STJ_19890312_XYZ.txt",
        ".scubas_config/1989/STJ_19890313_XYZ.txt",
        ".scubas_config/1989/STJ_19890314_XYZ.txt" 
    ]
)
stns, dates = (
    ["FRD", "STJ", "HAD"], 
    [dt.datetime(1989,3,13), dt.datetime(1989,3,14,12)]
)
frames = read_bfield_frames(stns, [FRD, HAD, STJ])

profiles = [
    PROFILES.CS_W, PROFILES.DO_1, PROFILES.DO_2, PROFILES.DO_3,
    PROFILES.DO_4, PROFILES.DO_5, PROFILES.MAR, PROFILES.DO_6,
    PROFILES.CS_E
]
tlines, cable = compile_cable_to_calculate_parameters(FRD, STJ, HAD, profiles)
omni = _load_omni_(dates)
ts = TimeSeriesPlot(
    [
        dt.datetime(1989,3,13), 
        dt.datetime(1989,3,14)
    ],
    fig_title="Date: 12-15 March 1989", text_size=12,
)
ts.add_omni(
    omni, ylelft="SymH", yright="AE", ylims=[[-600, 100], [0, 3000]],
    ytitles = [r"SymH [nT]", r"AE [nT]"],
)
ts.add_mag(frames)
ts.add_voltage(cable.tot_params)
ts.save("figures/Pot03.png")
ts.close()
# plot_potential(
#     frames, [dt.datetime(1989,3,12), dt.datetime(1989,3,15)],
#     cable.tot_params.index, [cable.tot_params["V(v)"]], omni,
#     major_locator=mdates.HourLocator(byhour=range(0, 24, 6)),
#     minor_locator=mdates.HourLocator(byhour=range(0, 24, 3)),
#     #yticks=[[-200, 0, 200, 400]], ylims=[[-200, 400], [-150, 150]],
#     formatter=DateFormatter("$%H^{%M}$"),
#     xlabel="Time [UT]", fig_title="Date: 12-15 March 1989", text_size=12,
#     fname="figures/Pot03.png"
# )

# plot_potential(
#     frames, [dt.datetime(2024,5,10,22), dt.datetime(2024,5,11)],
#     cable.tot_params.index, [cable.tot_params["V(v)"]], omni,
#     major_locator=mdates.HourLocator(byhour=range(0, 24, 1)),
#     minor_locator=mdates.MinuteLocator(byminute=range(0, 60, 15)),
#     #yticks=[[-200, 0, 200, 400]], ylims=[[-200, 400], [-150, 150]],
#     formatter=DateFormatter("%H^{%M}"),
#     xlabel="Minutes since 22 UT", fig_title="Date: 10 May 2024", text_size=12,
#     fname="figures/Pot02.png"
# )