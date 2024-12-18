import xarray as xr
import datetime as dt
from loguru import logger
import numpy as np
import aacgmv2
import pandas as pd
import glob
import os
import sys
sys.path.append("py/")
from plots import TimeSeriesPlot
from utils import get_cable_informations
from fetch_data import clean_B_fields

folder="simulation/May2024/"
base_folder="dataset/May2024/"
StationMap = dict(
    frd=dict(lat=38.3032, lon=-77.460),
    stj=dict(lat=47.56, lon=-52.71283),
    had=dict(lat=51, lon=-4.48),
)

def fetch_dataset(date, mlat, mlt):
    global folder, base_folder
    file = os.path.join(
        base_folder,
        f"{date.strftime('%Y%m%d')}.north.schavec-aacgm-supermag.60s.rev-0006.ncdf"
    )
    obj = xr.open_dataset(file, engine="netcdf4")
    dates = [
        dt.datetime(y,m,d,H,M,S)
        for y,m,d,H,M,S in zip(
            obj.variables["time_yr"].values.astype(int),
            obj.variables["time_mo"].values.astype(int),
            obj.variables["time_dy"].values.astype(int),
            obj.variables["time_hr"].values.astype(int),
            obj.variables["time_mt"].values.astype(int),
            obj.variables["time_sc"].values.astype(int),
        )
    ]
    i = dates.index(date)
    o = pd.DataFrame()
    o["mlat"], o["mlt"] = (
        obj.variables["mlat"].values[i, :],
        obj.variables["mlt"].values[i, :]
    )
    o["dbn_nez"] = obj.variables["dbn_nez"].values[i,:]
    o["dbe_nez"] = obj.variables["dbe_nez"].values[i,:]
    o["dbz_nez"] = obj.variables["dbz_nez"].values[i,:]
    o["db_nez"] = np.sqrt(o["dbn_nez"]**2+o["dbe_nez"]**2)
    mlt_hi, mlt_low = np.mod(mlt+1, 24), np.mod(mlt-1, 24)
    if mlt_hi < mlt_low:
        d = mlt_low
        mlt_low = mlt_hi
        mlt_hi = d
    x = o[
        (o.mlt>=mlt_low)
        & (o.mlt<=mlt_hi)
        & (o.mlat>=mlat-2)
        & (o.mlat<=mlat+2)
    ]
    x["dates"] = date
    logger.info(f"{date}, {np.max(x.db_nez)}")
    if len(x):
        x = x.groupby("dates").mean().reset_index()
    else:
        logger.info(f"No data point {date}/{mlat}/{mlt}")
        x["dates"] = [date]
        x["mlat"], x["mlt"] = [mlat], [mlt]
        x["dbn_nez"], x["dbe_nez"], x["dbz_nez"], x["db_nez"] = (
            [np.nan], [np.nan], [np.nan], [np.nan]
        )
    return x

def fetch_data_by_station(stn = "frd"):
    global folder, base_folder, StationMap
    os.makedirs(folder, exist_ok=True)
    location = StationMap[stn]
    frame = pd.DataFrame()
    fname = f"{folder}station_{stn}.csv"
    dates = [dt.datetime(2024,5,10) + dt.timedelta(minutes=t) for t in range(1440*2)]
    if not os.path.exists(fname):
        glat, glon = location["lat"], location["lon"]
        for date in dates:
            mlat, _, mlt = aacgmv2.get_aacgm_coord(glat, glon, 300, date)
            x = fetch_dataset(date, mlat, mlt)
            frame = pd.concat([frame, x])
        frame.to_csv(fname, float_format="%g", index=False, header=True)
    stn_validation_plots(stn)
    return

def stn_validation_plots(stn, dates=[dt.datetime(2024,5,10,12), dt.datetime(2024,5,12)]):
    global folder, base_folder
    from supermag import SuperMAGGetData, sm_grabme
    (_, data) = SuperMAGGetData("shibaji7", dates[0], 3600*36, "all,delta=start,baseline=yearly", stn.upper())
    if len(data) > 0:
        for comp in ["N", "E", "Z"]:
            for cord in ["nez", "geo"]:
                data[comp + "_" + cord] = sm_grabme(data, comp, cord)
        data.drop(["N", "E", "Z"], axis=1, inplace=True)
        data.tval = data.tval.apply(
            lambda x: dt.datetime.utcfromtimestamp(x)
        )
        print(data)
    fname, file = (
        f"{folder}station_{stn}.csv",
        f"{base_folder}{stn}20240510psec.sec.csv"
    )
    df = pd.read_csv(fname, parse_dates=["dates"])
    o = clean_B_fields([stn], [[file]])[stn]
    o.reset_index(inplace=True)
    o.X = o.X - np.nanmean(o.X.iloc[:60*10])
    o.Y = o.Y - np.nanmean(o.Y.iloc[:60*10])
    o.Z = o.Z - np.nanmean(o.Z.iloc[:60*10])
    import matplotlib.dates as mdates
    ts = TimeSeriesPlot(
        dates, f"", num_subplots=3, text_size=15,
        major_locator=mdates.HourLocator(byhour=range(0, 24, 4)),
        minor_locator=mdates.HourLocator(byhour=range(0, 24, 1)),
    )
    print(df.head())
    ax = ts.add_curve(df.dates, -df.dbn_nez,)
    ts.add_curve(o.Date, o.Y, ax = ax, color="b")
    ts.add_curve(data.tval, -data.N_nez, ax = ax, color="r", ylabel=r"$N_{nez}(Y)$ [nT]", ylim=[-1000, 1000])
    ax = ts.add_curve(df.dates, df.dbe_nez,)
    ts.add_curve(o.Date, o.X, ax = ax, color="b")
    ts.add_curve(data.tval, data.E_nez, ax = ax, color="r", ylabel=r"$E_{nez}(X)$ [nT]", ylim=[-1000, 1000])
    
    # ax = ts.add_curve(df.dates, df.dbn_nez, ylim=[-150, 200])
    # ts.add_curve(data.tval, data.N_nez, ax = ax, color="r", ylabel=r"$N_{nez}$ [nT]", ylim=[-500, 500])
    # ts.add_curve(o.Date, o.Y, ax = ax, color="b")
    # ax = ts.add_curve(df.dates, df.dbe_nez, ylim=[-150, 200])
    # ts.add_curve(data.tval, data.E_nez, ax = ax, color="r", ylabel=r"$E_{nez}$ [nT]", ylim=[-500, 500])
    # ts.add_curve(o.Date, o.X, ax = ax, color="b")
    ts.save(f"figures/SuperMAG.validation.{stn}.png")
    return

def interpolate_smag():
    global folder, base_folder
    os.makedirs(folder, exist_ok=True)
    cable = get_cable_informations()
    dates = [dt.datetime(2024,5,10) + dt.timedelta(minutes=t) for t in range(2*1440)]
    for i, seg in enumerate(cable.cable_seg):
        print(seg["center"]["lat"], seg["center"]["lon"])
        fname = f"{folder}cable_segment_{i}.csv"
        if not os.path.exists(fname):
            glat, glon = seg["center"]["lat"], seg["center"]["lon"]
            frame = pd.DataFrame()
            logger.info(f"Segment {i}")
            for date in dates:
                mlat, _, mlt = aacgmv2.get_aacgm_coord(glat, glon, 300, date)
                x = fetch_dataset(date, mlat, mlt)
                frame = pd.concat([frame, x])
            frame.rename(columns=dict(dates="Date",dbn_nez="X",dbe_nez="Y",dbz_nez="Z",db_nez="F"), inplace=True)
            frame.to_csv(fname, float_format="%g", index=False, header=True)
    return

def create_smag_stack_plot(dates=[dt.datetime(2024,5,10,12), dt.datetime(2024,5,12)]):
    global folder, base_folder
    import matplotlib.dates as mdates
    ts = TimeSeriesPlot(
        dates, f"", num_subplots=3, text_size=15,
        major_locator=mdates.HourLocator(byhour=range(0, 24, 4)),
        minor_locator=mdates.HourLocator(byhour=range(0, 24, 1)),
    )
    files = glob.glob(f"{folder}/cable_segment*.csv")
    files.sort()
    od = {}
    for i, f in enumerate(files):
        d = pd.read_csv(f, parse_dates=["dates"])
        d.set_index("dates", inplace=True)
        d = d.interpolate(method="linear")
        d.reset_index(inplace=True)
        d = d.rename(columns=dict(dbn_nez="Y", dbe_nez="X", dbz_nez="Z"))
        od[f"CS_{i+1}"] = d
        logger.info(f"File {f}/{len(d)}")
        print(d.head())
    ts.add_mag(od, stations=["CS_1"])
    ts.add_mag(od, stations=["CS_4"])
    ts.add_mag(od, stations=["CS_9"])
    # ts.add_mag(od, stations=["CS_0", "CS_1", "CS_2"])
    # ts.add_mag(od, stations=["CS_3", "CS_4", "CS_5"])
    # ts.add_mag(od, stations=["CS_6", "CS_7", "CS_8"])
    ts.save("figures/SuperMAG.stack.png")
    return

from cable import SCUBASModel

def run_May2024_storm():
    dates = [dt.datetime(2024,5,10,12), dt.datetime(2024,5,12)]
    segment_files = [
        ["simulation/May2024/cable_segment_0.csv"],
        ["simulation/May2024/cable_segment_1.csv"],
        ["simulation/May2024/cable_segment_2.csv"],
        ["simulation/May2024/cable_segment_3.csv"],
        ["simulation/May2024/cable_segment_4.csv"],
        ["simulation/May2024/cable_segment_5.csv"],
        ["simulation/May2024/cable_segment_6.csv"], 
        ["simulation/May2024/cable_segment_7.csv"],
        ["simulation/May2024/cable_segment_8.csv"],
    ]
    model = SCUBASModel(segment_files=segment_files)
    model.initialize_TL()
    model.run_cable_segment()
    data = model.cable.tot_params.reset_index().copy()
    data.to_csv("simulation/May2024/SCUBAS-Simulation-SuperMAG-Fitted.csv", float_format="%g", index=False)
    return

def SuperMAG_compare_plots():
    smag_simulations = pd.read_csv("simulation/May2024/SCUBAS-Simulation-SuperMAG-Fitted.csv", parse_dates=["Time"]).set_index("Time")
    print(smag_simulations.head())
    scuba_simulations = pd.read_csv("simulation/May2024/SCUBAS-Simulation-3-Stations.csv", parse_dates=["Time"]).set_index("Time")
    dates = [dt.datetime(2024,5,10,12), dt.datetime(2024,5,12)]
    import matplotlib.dates as mdates
    ts = TimeSeriesPlot(
        dates,
        major_locator=mdates.HourLocator(byhour=range(0, 24, 6)),
        minor_locator=mdates.HourLocator(byhour=range(0, 24, 3)),
        fig_title="Compare Voltages (Red: 3 Stations, Black: SuperMAG Fitted)", 
        text_size=15,
        num_subplots=1,
    )
    ax = ts.add_voltage(smag_simulations, xlabel="")
    ts.add_voltage(scuba_simulations, ax=ax, color="r", xlabel="Hours Since 10 May 2024, 12 UT")
    ts.save("simulation/May2024/SCUBAS-Simulation-SMag-Compare.png")
    ts.close()
    return

if __name__ == "__main__":
    # 1. Fetch o, dates based on 10 and 11 th May (D)
    # 2. Call fetch_dataset_from_netcdf_by_date_range for whole date range 
    #       and save to local files. (D)
    # TODO: 3. May be plan to add 12th May
    # interpolate_smag()
    fetch_data_by_station("frd")
    fetch_data_by_station("had")
    fetch_data_by_station("stj")

    # TODO: Plot all the supermag dataset for 9 different segments TS plot
    # create_smag_stack_plot()

    # TODO: Invoke scubas by the interpolated dataset 
    # run_May2024_storm()
    # SuperMAG_compare_plots()