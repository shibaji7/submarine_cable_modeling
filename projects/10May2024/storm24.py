from scubas.conductivity import ConductivityProfile
from scubas.datasets import Site, PROFILES
from scubas.models import OceanModel
from scubas.plotlib import plot_transfer_function, potential_along_section, cable_potential, update_rc_params
from scubas.cables import TransmissionLine, Cable
from scubas.conductivity import ConductivityProfile as CP

import pandas as pd
import datetime as dt
import numpy as np



def compile_cable_to_calculate_parameters(FRD_files, STJ_files, HAD_files, profiles):
    land50 = PROFILES.CS_E
    land50.layers[0].thickness = 50
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
                site=profiles[0],
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
                site=profiles[1],
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
                site=profiles[2],
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
                site=profiles[3],
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
                site=profiles[4],
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
                site=profiles[5],
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
                site=profiles[6],
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
                site=profiles[7],
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
                site=profiles[8],
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

def read_iaga(file, return_xyzf=True, return_header=False):
    """
    Read IAGA profiles
    """

    # Read Headers
    header_records = {"header_length": 0}

    with open(file, "r") as openfile:
        for newline in openfile:
            if newline[0] == " ":
                header_records["header_length"] += 1
                label = newline[1:24].strip()
                description = newline[24:-2].strip()
                header_records[label.lower()] = description

    if len(header_records["reported"]) % 4 != 0:
        raise ValueError(
            "The header record does not contain 4 values: {0}".format(
                header_records["reported"]
            )
        )
    record_length = len(header_records["reported"]) // 4
    column_names = [
        x for x in header_records["reported"][record_length - 1 :: record_length]
    ]
    seen_count = {}
    for i, col in enumerate(column_names):
        if col in seen_count:
            column_names[i] += str(seen_count[col])
            seen_count[col] += 1
        else:
            seen_count[col] = 1
    df = pd.read_csv(
        file,
        header=header_records["header_length"],
        delim_whitespace=True,
        parse_dates=[[0, 1]],
        infer_datetime_format=True,
        index_col=0,
        usecols=[0, 1, 3, 4, 5, 6],
        na_values=[99999.90, 99999.0, 88888.80, 88888.00],
        names=["Date", "Time"] + column_names,
    )
    df.index.name = "Date"
    if return_xyzf and "X" not in column_names and "Y" not in column_names:
        # Convert the data to XYZF format
        # Only convert HD
        if "H" not in column_names or "D" not in column_names:
            raise ValueError(
                "Only have a converter for HDZF->XYZF\n"
                + "Input file is: "
                + header_records["reported"]
            )

        # IAGA-2002 D is reported in minutes of arc.
        df["X"] = df["H"] * np.cos(np.deg2rad(df["D"] / 60.0))
        df["Y"] = df["H"] * np.sin(np.deg2rad(df["D"] / 60.0))
        del df["H"], df["D"]
    if return_header:
        return df, header_records
    else:
        return df

def read_Bfield_data(files, return_xyzf=True, csv_file_date_name="Date"):
    """
    Read B-Files
    """
    Bfield = pd.DataFrame()
    for file in files:
        file_type = file.split(".")[-1]
        if file_type == "txt":
            o = read_iaga(file, return_xyzf, return_header=False)
        elif file_type == "csv":
            o = pd.read_csv(file, parse_dates=[csv_file_date_name])
            o = o.rename(columns={csv_file_date_name: "Date"})
            o = o.set_index("Date")
            o.index.name = "Date"
        Bfield = pd.concat([Bfield, o])
    return Bfield

FRD, HAD, STJ = (
    [
        "mag/frd20240510psec.sec.txt"
    ], [
        "mag/had20240510psec.sec.txt"
    ], [
        "mag/stj20240510vsec.sec.txt"
    ]
)
import os
os.system("rm -rf mag/*")
os.system("cp bkp/* mag/")

stns = ["FRD", "STJ", "HAD"]
frames = dict()
for stn, fs in zip(stns, [FRD, STJ, HAD]):
    o = pd.DataFrame()
    o = pd.concat([o, read_Bfield_data(fs)])
    # Remove Datagaps
    print("Pre-Is there any issue / nan data? (X,Y,Z)", o.X.hasnans, o.Y.hasnans, o.Z.hasnans)
    o = o.replace(99999.0, np.nan)
    for key in ["X", "Y", "Z"]:                    
        o[key] = o[key].interpolate(method="from_derivatives")
        o[key] = o[key].ffill()
        o[key] = o[key].bfill()
    print("Post-Is there any issue / nan data? (X,Y,Z)", o.X.hasnans, o.Y.hasnans, o.Z.hasnans)
    fs[0] = fs[0].replace(".txt", ".csv")
    o.to_csv(fs[0], header=True, index=True, float_format="%g")
    frames[stn] = o

profiles = [
    PROFILES.CS_W, PROFILES.DO_1, PROFILES.DO_2, PROFILES.DO_3,
    PROFILES.DO_4, PROFILES.DO_5, PROFILES.MAR, PROFILES.DO_6,
    PROFILES.CS_E
]
tlines, cable = compile_cable_to_calculate_parameters(FRD, STJ, HAD, profiles)

import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
size = 12
mpl.rcParams.update(
        {"xtick.labelsize": size, "ytick.labelsize": size, "font.size": size}
    )


fig, axes = plt.subplots(nrows=2, ncols=1, dpi=300, figsize=(6, 4), sharex=True)

ax = axes[0]
ax.xaxis.set_major_formatter(DateFormatter(r"%H UT"))
ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 1)))
ax.set_ylabel(r"$B_{x,y}$ [nT]")
ax.set_ylim(-1800, 900)
for c, stn in zip(["k", "b", "g"], ["FRD", "STJ", "HAD"]):
    o = frames[stn]
    o.X = o.X - np.nanmean(o.X.iloc[:60*10])
    o.Y = o.Y - np.nanmean(o.Y.iloc[:60*10])
    o.Z = o.Z - np.nanmean(o.Z.iloc[:60*10])
    ax.plot(
        o.index, o.X, 
        color=c, ls="-", 
        lw=0.9, alpha=0.7,
        label=fr"$B[{stn}]$"
    )
    ax.plot(
        o.index, o.Y, 
        color=c, ls="--", 
        lw=0.9, alpha=0.7,
    )
ax.set_xlim(dt.datetime(2024,5,10,12), dt.datetime(2024,5,12))
ax.legend(loc=4)
ax.text(0.05, 1.05, "Date: 10-12 May 2024", ha="left", va="bottom", transform=ax.transAxes)


ax = axes[1]
ax.set_xlim(dt.datetime(2024,5,10,12), dt.datetime(2024,5,12))
ax.xaxis.set_major_formatter(DateFormatter(r"%H UT"))
ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 1)))
ax.set_ylabel("Potentials [Volts]")
ax.set_ylim(-500, 500)
ax.plot(
    cable.tot_params.index, 
    cable.tot_params["V(v)"], 
    color="r",
    ls="-", lw=1.5, 
    label=r"$\mathcal{E}_C$", 
    alpha=0.7
)
ax.set_xlim(dt.datetime(2024,5,10,12), dt.datetime(2024,5,12))
ax.legend(loc=1)
ax.set_xlabel("Time [UT]")
fig.savefig("figures/Pot.png", bbox_inches="tight")