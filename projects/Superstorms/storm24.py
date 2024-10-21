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
    _load_omni_,
    load_speadas
)
from event_plots import TimeSeriesPlot


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

# ['thc_fgs_btotal', 
# 'thc_fgs_gse', 'thc_fgs_gsm', 
# 'thc_fgs_dsl', 'thc_fgl_btotal', 'thc_fgl_gse', 
# 'thc_fgl_gsm', 'thc_fgl_dsl', 'thc_fgl_ssl', 'thc_fgh_btotal', 
# 'thc_fgh_gse', 'thc_fgh_gsm', 'thc_fgh_dsl', 'thc_fgh_ssl']
themis_fgm, themis_mom = load_speadas(dates)

tlines, cable = compile_cable_to_calculate_parameters(FRD, STJ, HAD, profiles)

ts = TimeSeriesPlot(
    [
        dt.datetime(2024,5,10,16,50), 
        dt.datetime(2024,5,10,17,20)
    ],
    major_locator=mdates.MinuteLocator(byminute=range(0, 60, 5)),
    minor_locator=mdates.MinuteLocator(byminute=range(0, 60, 1)),
    fig_title="Date: 10 May 2024 / IP Shock", text_size=15,
    num_subplots=3,
)
ts.add_vlines(
    ts.add_themis(themis_fgm, themis_mom, ["thc_fgs_gsm", "pdyn"]),
    vlines=[dt.datetime(2024,5,10,17,0,30)], colors=["r"]
)
ts.add_vlines(
    ts.add_mag(frames, ylim=[-200, 500]), 
    vlines=[dt.datetime(2024,5,10,17,6)], colors=["r"]
)
ts.add_vlines(
    ts.add_voltage(cable.tot_params, xlabel="Minutes since 17 UT", ylim=[-200, 100]),
    vlines=[dt.datetime(2024,5,10,17,6)], colors=["r"]
)
ts.save("figures/Pot01.png")
ts.close()

ts = TimeSeriesPlot(
    [
        dt.datetime(2024,5,10,21,30), 
        dt.datetime(2024,5,11)
    ],
    major_locator=mdates.HourLocator(byhour=range(0, 24, 1)),
    minor_locator=mdates.MinuteLocator(byminute=range(0, 60, 15)),
    fig_title="Date: 10-11 May 2024 / IMF By turning", text_size=15,
    num_subplots=3,
)
ts.add_vlines(
    ts.add_themis(themis_fgm, themis_mom, ["thc_fgs_gsm", "pdyn"]),
    vlines=[dt.datetime(2024,5,10,22,30)], colors=["k"]
)
ts.add_vlines(
    ts.add_mag(frames, ylim=[-1000, 1000]),
    vlines=[dt.datetime(2024,5,10,22,36)], colors=["k"]
)
ts.add_vlines(
    ts.add_voltage(cable.tot_params, xlabel="Minutes since 22 UT", ylim=[-50, 400]),
    vlines=[dt.datetime(2024,5,10,22,36)], colors=["k"]
)
ts.save("figures/Pot02.png")
ts.close()


    