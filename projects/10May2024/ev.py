import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
from scubas.datasets import Site, PROFILES
from scubas.models import OceanModel
from scubas.plotlib import plot_transfer_function, potential_along_section, cable_potential, update_rc_params
from scubas.cables import TransmissionLine, Cable
from scubas.conductivity import ConductivityProfile as CP

import pandas as pd
import datetime as dt
import numpy as np

import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

def compile_cable_to_calculate_parameters(FRD_files, STJ_files, HAD_files):
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
    print(fs)
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
print(fs)

mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
fig, axes = plt.subplots(nrows=2, ncols=1, dpi=240, figsize=(6, 6))
multiplier, colors = [1, 0, -1], ["r", "k", "b"]
base=1000

for i, stn in enumerate(["FRD", "STJ", "HAD"]):
    frame = frames[stn]
    print(stn, frame.head())
    ax = axes[0]
    ax.set_xlim(dt.datetime(2024,5,10), dt.datetime(2024,5,12))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
    ax.plot(frame.index, frame.X-np.mean(frame.X), colors[i], ls="-", lw=1., label=stn.upper())
    ax.set_ylabel(r"$B_x$, nT", fontdict={"color": "k"})
    ax.set_ylim(-1500,1000)
    #ax.axvline(frame.index.tolist()[1000], ymin=4/6, ymax=5/6, color = "darkgreen", drawstyle="steps-mid")
    #ax.axhline(2000, xmin=0.86, xmax=0.86+2e-2, color = "darkgreen")
    #ax.axhline(1000, xmin=0.86, xmax=0.86+2e-2, color = "darkgreen")
    #ax.text(frame.index.tolist()[-970], 1500, "1000 nT", ha="left", va="center", fontdict={"color": "darkgreen", "size":10})
    #ax.set_yticklabels([])
    ax.legend(loc=2)
    ax.text(0.95,0.95,"(a)",ha="right",va="center",transform=ax.transAxes)
    ax = axes[1]
    ax.set_xlim(dt.datetime(2024,5,10), dt.datetime(2024,5,12))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
    ax.plot(frame.index, frame.Y-np.mean(frame.Y), colors[i], ls="-", lw=1., label=stn.upper())
    ax.set_ylabel(r"$B_y$, nT", fontdict={"color": "k"})
    #ax.axvline(frame.index.tolist()[2000], ymin=4/6, ymax=5/6, color = "darkgreen", drawstyle="steps-mid")
    #ax.axhline(2000, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
    #ax.axhline(1000, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
    #ax.text(frame.index.tolist()[1970], 1500, "1000 nT", ha="right", va="center", fontdict={"color": "darkgreen", "size":10})
    ax.set_ylim(-500,500)
    #ax.set_xticklabels([])
    ax.text(0.95,0.95,"(b)",ha="right",va="center",transform=ax.transAxes)
axes[1].set_xlabel("Time, UT")
fig.subplots_adjust(wspace=.2, hspace=.1)
fig.savefig("figures/Bxy.Field.png", bbox_inches="tight")

tlines, cable = compile_cable_to_calculate_parameters(FRD, STJ, HAD)
o = cable.tot_params.copy().reset_index()
o.to_csv("elc/10May2024.csv", header=True, index=False, float_format="%g")

params = cable.tot_params
segment_name = ["CS-W", "DO-1", "DO-2", "DO-3", "DO-4", "DO-5", "MAR", "DO-6", "CS-E"]
stations = ["FRD", "FRD", "STJ", "STJ", "STJ", "STJ", "STJ", "HAD", "HAD"]

fig, axes = plt.subplots(nrows=2, ncols=1, dpi=240, figsize=(5, 8), sharex=True)
multiplier, colors = np.array([4,3,2,1,0,-1,-2,-3,-4])*3.5, ["r", "k", "b"]
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
          "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
base=90
for i in range(9):
    ax = axes[0]
    ax.text(0.05, 1.05, "Date: 10-12 May 2024", ha="left", va="bottom", transform=axes[0].transAxes)
    ax.set_xlim(dt.datetime(2024,5,10), dt.datetime(2024,5,12))
    ax.xaxis.set_major_formatter(DateFormatter(r"%H UT"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
    # ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    # ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
    ax.plot(
        params.index, 
        (base*multiplier[i])+params["E.X.%02d"%(i)]-np.mean(params["E.X.%02d"%(i)]), 
        colors[i], 
        ls="-", 
        lw=0.6
    )
    ax.set_xticklabels([])
    ax.set_ylabel(r"$E_x$ [mV/km]", fontdict={"color": "k"})
    ax.set_ylim(-1500,1500)
    ax.text(-0.06,0.2,"(a)",ha="left",va="bottom",transform=ax.transAxes)
    # if i==0:
    #     ax.axvline(params.index.tolist()[-1100], ymin=18.5/30, ymax=20.5/30, 
    #                color = "darkgreen", drawstyle="steps-mid")
    #     ax.axhline(350, xmin=0.82, xmax=0.82+1e-2, color = "darkgreen")
    #     ax.axhline(550, xmin=0.82, xmax=0.82+1e-2, color = "darkgreen")
    #     ax.text(params.index.tolist()[2800], 425, "200 mV/km", ha="left", va="center", 
    #             fontdict={"color": "darkgreen", "size":10})
    ax.set_yticklabels([])
    ax.legend(loc=2)
    ax = axes[1]
    ax.set_xlim(dt.datetime(2024,5,10), dt.datetime(2024,5,12))
    ax.xaxis.set_major_formatter(DateFormatter(r"%H UT"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
    # ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    # ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
    txt = r"%s[%s]"%(segment_name[i], stations[i])
    ax.plot(
        params.index, 
        (base*multiplier[i])+params["E.Y.%02d"%(i)]-np.mean(params["E.Y.%02d"%(i)]), 
        colors[i], 
        ls="-", 
        lw=0.6, 
        label=txt
    )
    ax.set_ylabel(r"$E_y$ [mV/km]", fontdict={"color": "k"})
    # if i==0:
    #     ax.axvline(params.index.tolist()[-1100], ymin=18.5/30, ymax=20.5/30, 
    #                color = "darkgreen", drawstyle="steps-mid")
    #     ax.axhline(350, xmin=0.82, xmax=0.82+1e-2, color = "darkgreen")
    #     ax.axhline(550, xmin=0.82, xmax=0.82+1e-2, color = "darkgreen")
    #     ax.text(params.index.tolist()[2800], 425, "200 mV/km", ha="left", va="center", 
    #             fontdict={"color": "darkgreen", "size":10})
        #ax.axvline(frame.index.tolist()[2000], ymin=4/6, ymax=5/6, color = "darkgreen", drawstyle="steps-mid")
        #ax.axhline(2000, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
        #ax.axhline(1000, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
        #ax.text(frame.index.tolist()[1970], 1500, "1000 nT", ha="right", va="center", fontdict={"color": "darkgreen", "size":10})
    #ax.set_ylim(-1500,1500)
    ax.set_yticklabels([])
    ax.text(-0.06,0.2,"(b)",ha="left",va="bottom",transform=ax.transAxes)
axes[1].set_xlabel("Time [UT]")
axes[1].legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=10)
fig.subplots_adjust(wspace=.1, hspace=.1)
fig.savefig("figures/EField.png", bbox_inches="tight")


fig, axes = plt.subplots(nrows=4, ncols=1, dpi=240, figsize=(6, 8), sharex=True)
for ax in axes:
    ax.set_xlim(dt.datetime(2024,5,10), dt.datetime(2024,5,12))
    ax.xaxis.set_major_formatter(DateFormatter(r"%H UT"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 1)))
    ax.set_ylabel("Potentials [Volts]")
    ax.set_ylim(-500, 500)
axes[0].plot(params.index, params["V(v)"], "darkgreen", ls="-", lw=1.5, label=r"$\mathcal{E}_C$", alpha=0.7)
axes[0].legend(loc=1)
axes[1].plot(params.index, params["U0"], "r", ls="-", lw=2, alpha=0.4, label=r"$U_W$")
axes[1].legend(loc=1)
axes[2].plot(params.index, params["U1"], "b", ls="-", lw=2, alpha=0.4, label=r"$U_E$")
axes[2].legend(loc=1)
axes[3].plot(params.index, params["V(v)"]+(params["U0"]-params["U1"]), "k", ls="-", lw=1., label=r"$V_{TAT-8}$")
axes[3].legend(loc=1)
axes[3].set_xlabel("Time [UT]")
axes[0].text(0.05, 1.05, "Date: 10-12 May 2024", ha="left", va="bottom", transform=axes[0].transAxes)
fig.savefig("figures/Pot.png", bbox_inches="tight")
#ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
