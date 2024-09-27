from scubas.conductivity import ConductivityProfile as CP
from event_runs import (
    read_bfield_frames, initialize_mag_fold, 
    run_mc, get_dates, read_mc_datasets, 
    get_mcmc_outputs_CI,
    extract_mc_transfer_functions_by_section
)
from event_plots import (
    plot_potential,
    plot_transfer_functions
)

import matplotlib.dates as mdates
import datetime as dt
import numpy as np

initialize_mag_fold()
FRD, HAD, STJ = (
    [".scubas_config/frd20240510psec.sec.txt"],
    [".scubas_config/had20240510psec.sec.txt"], 
    [".scubas_config/stj20240510vsec.sec.txt"]
)
stns, dates, N = (
    ["FRD", "STJ", "HAD"], 
    [dt.datetime(2024,5,10,17), dt.datetime(2024,5,10,18)],
    10
)
frames = read_bfield_frames(stns, [FRD, HAD, STJ])

latlons = [
    dict(lat=39.6, lon=-74.33), dict(lat=38.79, lon=-72.62),
    dict(lat=37.11, lon=-68.94), dict(lat=39.80, lon=-48.20),
    dict(lat=40.81, lon=-45.19), dict(lat=43.15, lon=-39.16),
    dict(lat=44.83, lon=-34.48), dict(lat=46.51, lon=-22.43),
    dict(lat=47.85, lon=-9.05), dict(lat=50.79, lon=-4.55),
]
binned_lat_lon = []
for l in range(len(latlons)):   
    binned_lat_lon.append(
        [latlons[l]["lat"], latlons[l]["lon"]]
    )
mc_profiles = CP.compile_bined_profiles(np.array(binned_lat_lon), n=N)


keys_to_store = [
    "V(v)", "Vt(v)", 
    # 'E.Y', 'E.Y.00', 'E.Y.01', 'E.Y.02', 'E.Y.03', 'E.Y.04', 'E.Y.05', 'E.Y.06',
    # 'E.Y.07', 'E.Y.08', 'E.X', 'E.X.00', 'E.X.01', 'E.X.02', 'E.X.03', 'E.X.04',
    # 'E.X.05', 'E.X.06', 'E.X.07', 'E.X.08', 'V(v).00', 'V(v).01', 'V(v).02', 'V(v).03',
    # 'V(v).04', 'V(v).05', 'V(v).06', 'V(v).07', 'V(v).08', 'U0', 'U1'
]

run_mc(
    mc_profiles, keys_to_store,
    FRD, STJ, HAD, n=N
)
df = get_mcmc_outputs_CI(
    read_mc_datasets(n=N)
)

plot_potential(
    frames, dates,
    get_dates(
        [dt.datetime(2024,5,10), dt.datetime(2024,5,12)]
    ), 
    [
        df["mean"], 
        df["CI"][0]["ub"],
        df["CI"][0]["lb"],
        df["CI"][1]["ub"],
        df["CI"][1]["lb"],
    ],
    linewdiths=[0.7, 0.5, 0.4, 0.4],
    colors=[["k", "b", "g"], ["k", "b", "b", "r", "r"]],
    fname="figures/Pot_mc_runs.png",
    fig_title="Date: 17-18 UT, 10 May 2024",
    yticks=[[-500, -100, 0, 300]],
    ylims=[[-500, 300], [-300, 100]],
    major_locator=mdates.MinuteLocator(byminute=range(0, 60, 10)),
    minor_locator=mdates.MinuteLocator(byminute=range(0, 60, 10)),
    alphas=[1., 0.5, 0.5, 0.3, 0.3]
)

flim, M = [1e-6, 1e0], 100
freqs = np.linspace(
    flim[0], flim[1], 
    int(flim[1] / (M*flim[0])) + 1
)
df = extract_mc_transfer_functions_by_section(
    freqs, mc_profiles, n=N
)
plot_transfer_functions(
    df, yticks=[[1e-3, 1e-0, 1e3],[0, 15, 30, 45, 60]],
    ylims=[[1e-3,1e3],[0, 60]],
)