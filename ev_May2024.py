import matplotlib as mpl
import matplotlib.pyplot as plt

plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]
import datetime as dt

import matplotlib.dates as mdates
import numpy as np
import pandas as pd
from matplotlib.dates import DateFormatter


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
    df.index.name = "Time"
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
        print(file, file_type)
        if file_type == "txt":
            o = read_iaga(file, return_xyzf, return_header=False)
        elif file_type == "csv":
            o = pd.read_csv(file, parse_dates=[csv_file_date_name])
            o = o.rename(columns={csv_file_date_name: "Time"})
            o = o.set_index("Time")
            o.index.name = "Time"
        Bfield = pd.concat([Bfield, o])
    return Bfield


FRD, HAD, STJ = (
    ["input/data/2024/frd20240510psec.sec.txt"],
    ["input/data/2024/had20240510vsec.sec.txt"],
    ["input/data/2024/stj20240510vsec.sec.txt"],
)

stns = ["FRD", "STJ", "HAD"]
frames = dict()
for stn, fs in zip(stns, [FRD, STJ, HAD]):
    o = pd.DataFrame()
    o = pd.concat([o, read_Bfield_data(fs)])
    frames[stn] = o
    print(o.isnull().values.any())

mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize": 12, "font.size": 12})
fig, axes = plt.subplots(nrows=2, ncols=1, dpi=240, figsize=(6, 6))
multiplier, colors = [1, 0, -1], ["r", "k", "b"]
base = 1000

for i, stn in enumerate(["FRD", "STJ", "HAD"]):
    frame = frames[stn]
    ax = axes[0]
    ax.set_xlim(dt.datetime(2024, 5, 10), dt.datetime(2024, 5, 13))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
    ax.plot(
        frame.index,
        (base * multiplier[i]) + frame.X - np.mean(frame.X),
        colors[i],
        ls="-",
        lw=1.0,
        label=stn.upper(),
    )
    ax.set_ylabel(r"$B_x$, nT", fontdict={"color": "k"})
    ax.set_ylim(-1500, 1500)
    # ax.axvline(frame.index.tolist()[1000], ymin=4/6, ymax=5/6, color = "darkgreen", drawstyle="steps-mid")
    # ax.axhline(2000, xmin=0.86, xmax=0.86+2e-2, color = "darkgreen")
    # ax.axhline(1000, xmin=0.86, xmax=0.86+2e-2, color = "darkgreen")
    # ax.text(frame.index.tolist()[-970], 1500, "1000 nT", ha="left", va="center", fontdict={"color": "darkgreen", "size":10})
    # ax.set_yticklabels([])
    ax.legend(loc=2)
    ax.text(0.95, 0.95, "(a)", ha="right", va="center", transform=ax.transAxes)
    ax = axes[1]
    ax.set_xlim(dt.datetime(2024, 5, 10), dt.datetime(2024, 5, 13))
    ax.xaxis.set_major_formatter(DateFormatter("%b.%d"))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_formatter(DateFormatter("%H UT"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 12)))
    ax.plot(
        frame.index,
        (base * multiplier[i]) + frame.Y - np.mean(frame.Y),
        colors[i],
        ls="-",
        lw=1.0,
        label=stn.upper(),
    )
    ax.set_ylabel(r"$B_y$, nT", fontdict={"color": "k"})
    # ax.axvline(frame.index.tolist()[2000], ymin=4/6, ymax=5/6, color = "darkgreen", drawstyle="steps-mid")
    # ax.axhline(2000, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
    # ax.axhline(1000, xmin=0.25, xmax=0.25+2e-2, color = "darkgreen")
    # ax.text(frame.index.tolist()[1970], 1500, "1000 nT", ha="right", va="center", fontdict={"color": "darkgreen", "size":10})
    ax.set_ylim(-1500, 1500)
    # ax.set_xticklabels([])
    ax.text(0.95, 0.95, "(b)", ha="right", va="center", transform=ax.transAxes)
axes[1].set_xlabel("Time, UT")
fig.subplots_adjust(wspace=0.2, hspace=0.1)
fig.savefig("Bxy.Field.png", bbox_inches="tight")

# land50 = PROFILES.CS_E
# land50.layers[0].thickness = 50
# # Create cable by each cable sections
# tlines = []
# tlines.append(
#     TransmissionLine(
#         sec_id="CS-W",
#         directed_length=dict(
#             edge_locations=dict(
#                 initial=dict(lat=39.6, lon=-74.33),
#                 final=dict(lat=38.79, lon=-72.62)
#             )
#         ),
#         elec_params=dict(
#             site=PROFILES.CS_W,
#             width=1.0,
#             flim=[1e-6, 1e0],
#         ),
#         active_termination=dict(
#             right=None,
#             left=PROFILES.LD,
#         ),
#     ).compile_oml(FRD),
# )
# tlines.append(
#     TransmissionLine(
#         sec_id="DO-1",
#         directed_length=dict(
#             edge_locations=dict(
#                 initial=dict(lat=38.79, lon=-72.62),
#                 final=dict(lat=37.11, lon=-68.94)
#             )
#         ),
#         elec_params=dict(
#             site=PROFILES.DO_1,
#             width=1.0,
#             flim=[1e-6, 1e0],
#         ),
#         active_termination=dict(
#             right=None,
#             left=None,
#         ),
#     ).compile_oml(FRD)
# )
# tlines.append(
#     TransmissionLine(
#         sec_id="DO-2",
#         directed_length=dict(
#             edge_locations=dict(
#                 initial=dict(lat=37.11, lon=-68.94),
#                 final=dict(lat=39.80, lon=-48.20)
#             )
#         ),
#         elec_params=dict(
#             site=PROFILES.DO_2,
#             width=1.0,
#             flim=[1e-6, 1e0],
#         ),
#         active_termination=dict(
#             right=None,
#             left=None,
#         ),
#     ).compile_oml(STJ)
# )
# tlines.append(
#     TransmissionLine(
#         sec_id="DO-3",
#         directed_length=dict(
#             edge_locations=dict(
#                 initial=dict(lat=39.80, lon=-48.20),
#                 final=dict(lat=40.81, lon=-45.19)
#             )
#         ),
#         elec_params=dict(
#             site=PROFILES.DO_3,
#             width=1.0,
#             flim=[1e-6, 1e0],
#         ),
#         active_termination=dict(
#             right=None,
#             left=None,
#         ),
#     ).compile_oml(STJ)
# )
# tlines.append(
#     TransmissionLine(
#         sec_id="DO-4",
#         directed_length=dict(
#             edge_locations=dict(
#                 initial=dict(lat=40.81, lon=-45.19),
#                 final=dict(lat=43.15, lon=-39.16)
#             )
#         ),
#         elec_params=dict(
#             site=PROFILES.DO_4,
#             width=1.0,
#             flim=[1e-6, 1e0],
#         ),
#         active_termination=dict(
#             right=None,
#             left=None,
#         ),
#     ).compile_oml(STJ)
# )
# tlines.append(
#     TransmissionLine(
#         sec_id="DO-5",
#         directed_length=dict(
#             edge_locations=dict(
#                 initial=dict(lat=43.15, lon=-39.16),
#                 final=dict(lat=44.83, lon=-34.48)
#             )
#         ),
#         elec_params=dict(
#             site=PROFILES.DO_5,
#             width=1.0,
#             flim=[1e-6, 1e0],
#         ),
#         active_termination=dict(
#             right=None,
#             left=None,
#         ),
#     ).compile_oml(STJ)
# )
# tlines.append(
#     TransmissionLine(
#         sec_id="MAR",
#         directed_length=dict(
#             edge_locations=dict(
#                 initial=dict(lat=44.83, lon=-34.48),
#                 final=dict(lat=46.51, lon=-22.43)
#             )
#         ),
#         elec_params=dict(
#             site=PROFILES.MAR,
#             width=1.0,
#             flim=[1e-6, 1e0],
#         ),
#         active_termination=dict(
#             right=None,
#             left=None,
#         ),
#     ).compile_oml(STJ)
# )
# tlines.append(
#     TransmissionLine(
#         sec_id="DO-6",
#         directed_length=dict(
#             edge_locations=dict(
#                 initial=dict(lat=46.51, lon=-22.43),
#                 final=dict(lat=47.85, lon=-9.05)
#             )
#         ),
#         elec_params=dict(
#             site=PROFILES.DO_6,
#             width=1.0,
#             flim=[1e-6, 1e0],
#         ),
#         active_termination=dict(
#             right=None,
#             left=None,
#         ),
#     ).compile_oml(STJ)
# )
# tlines.append(
#     TransmissionLine(
#         sec_id="CS-E",
#         directed_length=dict(
#             edge_locations=dict(
#                 initial=dict(lat=47.85, lon=-9.05),
#                 final=dict(lat=50.79, lon=-4.55)
#             )
#         ),
#         elec_params=dict(
#             site=PROFILES.CS_E,
#             width=1.0,
#             flim=[1e-6, 1e0],
#         ),
#         active_termination=dict(
#             right=land50,
#             left=None,
#         ),
#     ).compile_oml(STJ)
# )


# cable = Cable(tlines, tlines[0].components)
# cable.tot_params.head()
# cable.tot_params.to_csv("2024Storm.csv")

# params = cable.tot_params
# segment_name = ["CS-W", "DO-1", "DO-2", "DO-3", "DO-4", "DO-5", "MAR", "DO-6", "CS-E"]
# stations = ["FRD", "FRD", "STJ", "STJ", "STJ", "STJ", "STJ", "HAD", "HAD"]
