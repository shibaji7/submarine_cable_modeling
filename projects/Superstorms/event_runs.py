import pandas as pd
import numpy as np
import datetime as dt
import glob
from scipy import constants as C

from scubas.datasets import PROFILES
from scubas.cables import TransmissionLine, Cable
from scubas.conductivity import ConductivityProfile as CP
import os

def stack_profiles_to_csv(profiles, file=".scubas_config/sample.csv"):
    import os
    os.remove(file)
    open(file, "a").close()
    for o in profiles:
        o[0].to_csv(file, mode="a+", index=False, header=True)
    return

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

def read_bfield_frames(stns, stn_files):
    frames = dict()
    for stn, fs in zip(stns, stn_files):
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
    return frames

def initialize_mag_fold():
    import os
    os.system("cp -r dataset/* .scubas_config/")
    return

def run_mc(
    mc_profiles, 
    keys_to_store,
    FRD, STJ, HAD,
    fname=".scubas_config/{}.csv.bz2",
    n=1000,
):
    for i, mc_profile in enumerate(mc_profiles[:n]):
        fn = fname.format("%03d"%i)
        print(f"F_name: {fn}")
        if not os.path.exists(fn):
            tlines, cable = compile_cable_to_calculate_parameters(FRD, STJ, HAD, mc_profile)
            o = cable.tot_params[keys_to_store]
            o.to_csv(fname, index=False, float_format="%g", compression="bz2")
            del tlines, cable
    return

def get_dates(dates, resolution_sec=1):
    seconds = int((dates[1] - dates[0]).total_seconds()/resolution_sec)
    date_list = [dates[0]+dt.timedelta(seconds=s) for s in range(seconds)]
    return date_list

def read_mc_datasets(fstr = ".scubas_config/*.bz2", pname="V(v)", n=1000):
    files = glob.glob(fstr)
    files.sort()
    data = []
    for file in files[:n]:
        o = pd.read_csv(file, compression="bz2")
        data.append(o[pname].tolist())
    return np.array(data)

def get_mcmc_outputs_CI(data, ci_multipliers=[1, 1.96], dx=1):
    mean = data.mean(axis=0)[::dx]
    std = data.std(axis=0)[::dx]
    print(f"Shape: {mean.shape}/{std.shape}")
    o = dict(
        mean = mean, std = std,
        CI = [
            dict(
                ub = mean + ci*std, 
                lb = mean - ci*std
            )
            for ci in ci_multipliers
        ]
    )
    return o

def extract_mc_sites_by_section(mc_profiles, section=0, n=1000):
    sites = [profiles[section] for profiles in mc_profiles[:n]]
    return sites

def extract_mc_transfer_functions_by_section(freqs, mc_profiles, section=0, n=1000):
    sites = extract_mc_sites_by_section(mc_profiles, section, n)
    tfs = []
    for site in sites:
        tfs.append(calcTF(site, freqs))
    amp_val, phase_val = (
        np.array([t.amp.tolist() for t in tfs]).T,
        np.array([t.phase.tolist() for t in tfs]).T
    )
    stats = pd.DataFrame()
    (
        stats["freqs"], stats["amp"], 
        stats["amp_ub_1"], stats["amp_lb_1"],
        stats["amp_ub_2"], stats["amp_lb_2"], stats["phase"], 
        stats["phase_ub_1"], stats["phase_lb_1"],
        stats["phase_ub_2"], stats["phase_lb_2"]
    ) = (
        freqs, amp_val.mean(axis=1), 
        amp_val.mean(axis=1) + amp_val.std(axis=1),
        amp_val.mean(axis=1) - amp_val.std(axis=1),
        amp_val.mean(axis=1) + 1.96*amp_val.std(axis=1),
        amp_val.mean(axis=1) - 1.96*amp_val.std(axis=1),
        phase_val.mean(axis=1), 
        phase_val.mean(axis=1) + phase_val.std(axis=1),
        phase_val.mean(axis=1) - phase_val.std(axis=1),
        phase_val.mean(axis=1) + 1.96*phase_val.std(axis=1),
        phase_val.mean(axis=1) - 1.96*phase_val.std(axis=1)
    )
    o = dict(
        tfs=tfs,
        stats=stats,
    )
    return o

def calcZ(site, freqs):
    omega = 2 * C.pi * freqs
    sigma_s = 1 / site.layers[0].resistivity
    k2 = 1.0j * omega * C.mu_0 * sigma_s
    k = np.sqrt(k2)
    Zo = (1.0j * omega * C.mu_0 / k) / (C.mu_0 / 1.0e-3)
    return Zo

def calcTF(site, freqs):
    """
    Calculate the transfer functions.
    """
    Zo = calcZ(site, freqs)
    Zd = site.calcZ(freqs)[1]
    omega = 2 * C.pi * freqs
    sigma_s = 1 / site.layers[0].resistivity
    k2 = 1.0j * omega * C.mu_0 * sigma_s
    k = np.sqrt(k2)
    kd = k * site.layers[0].thickness

    func = Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Zo))
    tf = pd.DataFrame()
    tf["freq"], tf["E2B"], tf["amp"], tf["phase"] = (
        freqs, func, np.abs(func), 
        np.angle(func, deg=True)
    )
    return tf