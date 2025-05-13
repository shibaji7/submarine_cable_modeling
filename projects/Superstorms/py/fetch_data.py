import datetime as dt
import os

import numpy as np
import pandas as pd
import pyspedas
from loguru import logger

os.environ["OMNIDATA_PATH"] = "/home/shibaji/omni/"


def _load_omni_(dates, res=1):
    import pyomnidata

    logger.info(f"OMNIDATA_PATH: {os.environ['OMNIDATA_PATH']}")
    pyomnidata.UpdateLocalData()
    omni = pd.DataFrame(pyomnidata.GetOMNI(dates[0].year, Res=res))
    omni["time"] = omni.apply(
        lambda r: (
            dt.datetime(
                int(str(r.Date)[:4]),
                int(str(r.Date)[4:6]),
                int(str(r.Date)[6:].replace(".0", "")),
            )
            + dt.timedelta(hours=r.ut)
        ),
        axis=1,
    )
    omni = omni[(omni.time >= dates[0]) & (omni.time <= dates[1])]
    return omni


def load_speadas(dates, probe="c"):
    time_range = [
        dates[0].strftime("%Y-%m-%d/%H:%M"),
        dates[1].strftime("%Y-%m-%d/%H:%M"),
    ]
    data_fgm = pyspedas.themis.fgm(
        probe=probe, trange=time_range, time_clip=True, no_update=False, notplot=True
    )
    data_mom = pyspedas.themis.mom(
        probe=probe, trange=time_range, notplot=True, no_update=False, time_clip=True
    )
    pdyn = {
        "x": data_mom["thc_peem_density"]["x"],
        "y": data_mom["thc_peem_density"]["y"]
        * 1.67
        * (10 ** (-6))
        * 0.5
        * np.nansum(data_mom["thc_peim_velocity_gse"]["y"] ** 2, axis=1),
    }
    data_mom["pdyn"] = pdyn
    return data_fgm, data_mom


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


def clean_B_fields(stns, stn_files):
    frames = dict()
    for stn, fs in zip(stns, stn_files):
        o = pd.DataFrame()
        o = pd.concat([o, read_Bfield_data(fs)])
        # Remove Datagaps
        print(
            "Pre-Is there any issue / nan data? (X,Y,Z)",
            o.X.hasnans,
            o.Y.hasnans,
            o.Z.hasnans,
        )
        o = o.replace(99999.0, np.nan)
        for key in ["X", "Y", "Z"]:
            o[key] = o[key].interpolate(method="from_derivatives")
            o[key] = o[key].ffill()
            o[key] = o[key].bfill()
        print(
            "Post-Is there any issue / nan data? (X,Y,Z)",
            o.X.hasnans,
            o.Y.hasnans,
            o.Z.hasnans,
        )
        fs[0] = fs[0].replace(".txt", ".csv")
        o.to_csv(fs[0], header=True, index=True, float_format="%g")
        frames[stn] = o
    return frames
