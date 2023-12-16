# %%
"""
## This code is related to the Extereme Value Analysis
"""

# %%
from scubas.datasets import PROFILES
from scubas.models import OceanModel
from scubas.plotlib import plot_transfer_function, potential_along_section, cable_potential, update_rc_params
from scubas.cables import TransmissionLine, Cable
from scubas.conductivity import ConductivityProfile as CP

land50 = PROFILES.CS_E
land50.layers[0].thickness = 50

def compile_cable_to_calculate_parameters(FRD_files, STJ_files, HAD_files):
    # Create cable by each cable sections
    tlines = []
    tlines.append(
        TransmissionLine(
            sec_id="CS-W",
            directed_length=dict(
                edge_locations=dict(
                    initial=dict(lat=39.6, lon=-74.33), 
                    final=dict(lat=38.79, lon=-72.62)
                )
            ),
            elec_params=dict(
                site=PROFILES.CS_W,
                width=1.0,
                flim=[1e-6, 1e0],
            ),
            active_termination=dict(
                right=None,
                left=PROFILES.LD,
            ),
        ).compile_oml(FRD_files),
    )
    tlines.append(
        TransmissionLine(
            sec_id="DO-1",
            directed_length=dict(
                edge_locations=dict(
                    initial=dict(lat=38.79, lon=-72.62), 
                    final=dict(lat=37.11, lon=-68.94)
                )
            ),
            elec_params=dict(
                site=PROFILES.DO_1,
                width=1.0,
                flim=[1e-6, 1e0],
            ),
            active_termination=dict(
                right=None,
                left=None,
            ),
        ).compile_oml(FRD_files)
    )
    tlines.append(
        TransmissionLine(
            sec_id="DO-2",
            directed_length=dict(
                edge_locations=dict(
                    initial=dict(lat=37.11, lon=-68.94), 
                    final=dict(lat=39.80, lon=-48.20)
                )
            ),
            elec_params=dict(
                site=PROFILES.DO_2,
                width=1.0,
                flim=[1e-6, 1e0],
            ),
            active_termination=dict(
                right=None,
                left=None,
            ),
        ).compile_oml(STJ_files)
    )
    tlines.append(
        TransmissionLine(
            sec_id="DO-3",
            directed_length=dict(
                edge_locations=dict(
                    initial=dict(lat=39.80, lon=-48.20), 
                    final=dict(lat=40.81, lon=-45.19)
                )
            ),
            elec_params=dict(
                site=PROFILES.DO_3,
                width=1.0,
                flim=[1e-6, 1e0],
            ),
            active_termination=dict(
                right=None,
                left=None,
            ),
        ).compile_oml(STJ_files)
    )
    tlines.append(
        TransmissionLine(
            sec_id="DO-4",
            directed_length=dict(
                edge_locations=dict(
                    initial=dict(lat=40.81, lon=-45.19), 
                    final=dict(lat=43.15, lon=-39.16)
                )
            ),
            elec_params=dict(
                site=PROFILES.DO_4,
                width=1.0,
                flim=[1e-6, 1e0],
            ),
            active_termination=dict(
                right=None,
                left=None,
            ),
        ).compile_oml(STJ_files)
    )
    tlines.append(
        TransmissionLine(
            sec_id="DO-5",
            directed_length=dict(
                edge_locations=dict(
                    initial=dict(lat=43.15, lon=-39.16), 
                    final=dict(lat=44.83, lon=-34.48)
                )
            ),
            elec_params=dict(
                site=PROFILES.DO_5,
                width=1.0,
                flim=[1e-6, 1e0],
            ),
            active_termination=dict(
                right=None,
                left=None,
            ),
        ).compile_oml(STJ_files)
    )
    tlines.append(
        TransmissionLine(
            sec_id="MAR",
            directed_length=dict(
                edge_locations=dict(
                    initial=dict(lat=44.83, lon=-34.48), 
                    final=dict(lat=46.51, lon=-22.43)
                )
            ),
            elec_params=dict(
                site=PROFILES.MAR,
                width=1.0,
                flim=[1e-6, 1e0],
            ),
            active_termination=dict(
                right=None,
                left=None,
            ),
        ).compile_oml(STJ_files)
    )
    tlines.append(
        TransmissionLine(
            sec_id="DO-6",
            directed_length=dict(
                edge_locations=dict(
                    initial=dict(lat=46.51, lon=-22.43), 
                    final=dict(lat=47.85, lon=-9.05)
                )
            ),
            elec_params=dict(
                site=PROFILES.DO_6,
                width=1.0,
                flim=[1e-6, 1e0],
            ),
            active_termination=dict(
                right=None,
                left=None,
            ),
        ).compile_oml(HAD_files)
    )
    tlines.append(
        TransmissionLine(
            sec_id="CS-E",
            directed_length=dict(
                edge_locations=dict(
                    initial=dict(lat=47.85, lon=-9.05), 
                    final=dict(lat=50.79, lon=-4.55)
                )
            ),
            elec_params=dict(
                site=PROFILES.CS_E,
                width=1.0,
                flim=[1e-6, 1e0],
            ),
            active_termination=dict(
                right=land50,
                left=None,
            ),
        ).compile_oml(HAD_files)
    )
    
    # Running the cable operation
    cable = Cable(tlines, tlines[0].components)
    return tlines, cable

# %%
# Download the datasets, save to local, and validate datagaps

import numpy as np
import shutil
import requests
import os
import datetime as dt
import pandas as pd

clear = False
MAG_FOLDER, ELC_FOLDER = "ev/mag/", "ev/elc/"
if clear: 
    shutil.rmtree(MAG_FOLDER)
    shutil.rmtree(ELC_FOLDER)
os.makedirs(MAG_FOLDER, exist_ok=True)
os.makedirs(ELC_FOLDER, exist_ok=True)

BASE_URL = "https://imag-data.bgs.ac.uk/GIN_V1/hapi"

years, magnetometers, coord_keys = (
    range(1997,2021), ["stj", "frd", "had"], ["X", "Y", "Z"]
)

for mag in magnetometers:
    for year in years:
        dates = [
            dt.datetime(year, 1, 1),
            dt.datetime(year+1, 1, 1)
        ]
        expected_dp = (dates[1]-dates[0]).total_seconds() / 60.
        dates_str = [
            dates[0].strftime("%Y-%m-%dT%H:%M:%SZ"),
            dates[1].strftime("%Y-%m-%dT%H:%M:%SZ")
        ]
        link = f"/data?dataset={mag}/definitive/PT1M/native&parameters=Field_Vector&start={dates_str[0]}&stop={dates_str[1]}"
        url = BASE_URL + link
        print("Downloadable article:", url)
        fname = MAG_FOLDER + f"{mag.upper()}.{dates[0].year}.csv"
        # Download and save to file
        if not os.path.exists(fname):
            r = requests.get(url)
            print(f"Response code <{mag.upper()}>:", r.status_code)
            txt = "Date,X,Y,Z\n" + r.text
            with open(fname, "w") as f:
                f.write(txt)
            # Validate number of data points
            o = pd.read_csv(fname, parse_dates=["Date"])
            print(f"Expecetd / Actual # data points: {expected_dp} / {len(o)}")
            # Remove Datagaps
            o = o.replace(99999.0, np.nan)
            for key in coord_keys:                    
                o[key] = o[key].interpolate(method="from_derivatives")
                o[key] = o[key].ffill()
                o[key] = o[key].bfill()
            print("Is there any issue / nan data? ", o.isnull().values.any())
            # Baseline removal
            for key in coord_keys:
                fmed = np.nanmedian(o[key][:120])
                o[key] = o[key] - fmed
            o.to_csv(fname, header=True, index=False, float_format="%g")

for year in years:            
    FRD_files, STJ_files, HAD_files = (
        [MAG_FOLDER + f"FRD.{year}.csv"],
        [MAG_FOLDER + f"STJ.{year}.csv"],
        [MAG_FOLDER + f"HAD.{year}.csv"]
    )
    fname = ELC_FOLDER + f"{year}.csv"
    if not os.path.exists(fname):
        tlines, cable = compile_cable_to_calculate_parameters(FRD_files, STJ_files, HAD_files)
        o = cable.tot_params.copy().reset_index()
        o.to_csv(fname, header=True, float_format="%g")
    else:
        o = pd.read_csv(fname, parse_dates=["Time"])
    o["hour"] = o.Time.apply(lambda x: x.hour + (x.dayofyear-1)*24)
    parameters = [
        "hour", "Time", "Vt(v)", "V(v)",
        "V(v).00", "V(v).01", "V(v).02",
        "V(v).03", "V(v).04", "V(v).05",
        "V(v).06", "V(v).07", "V(v).08",
        "U0", "U1", "E.Y", "E.Y.00", 
        "E.Y.01", "E.Y.02", "E.Y.03",
        "E.Y.04", "E.Y.05", "E.Y.06",
        "E.Y.07", "E.Y.08", "E.X", "E.X.00",
        "E.X.01", "E.X.02", "E.X.03",
        "E.X.04", "E.X.05", "E.X.06", "E.X.07", "E.X.08"
    ]
    o = o[parameters]
    for p in parameters[2:]:
        o[p] = o[p].abs()
    o["V_W"] = o["V(v).00"] + o["V(v).01"] + o["U0"]
    o["V_E"] = o["V(v).07"] + o["V(v).08"] - o["U1"]
    o["V_M"] = o["V(v).02"] + o["V(v).03"] + o["V(v).04"] + o["V(v).05"] + o["V(v).06"]
    o["V_SW"] = o["V(v).00"] + o["U0"]
    o["V_SE"] = o["V(v).08"] - o["U1"]
    o["V_D"] = o["V(v).01"] + o["V(v).02"] + o["V(v).03"] + o["V(v).04"] + o["V(v).05"] + o["V(v).06"] + o["V(v).07"]
    parameters.extend(["V_W", "V_E", "V_M", "V_SW", "V_SE", "V_D"])
    o = o.groupby(by="hour").agg([np.nanmin, np.nanmax])
    df = pd.DataFrame()
    df["start_time"], df["end_time"] = o.Time["nanmin"], o.Time["nanmax"]
    for p in parameters[2:]:
        (df[f"{p.replace('(v)', '')}_min"], df[f"{p.replace('(v)', '')}_max"]) = (
            o[p]["nanmin"], o[p]["nanmax"]
        )
    df.to_csv(fname.replace(".csv", ".hr.csv"), header=True, index=False, float_format="%g")