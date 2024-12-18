
import pandas as pd
import datetime as dt
from loguru import logger
import glob
import sys
import numpy as np
sys.path.append("py/")
from plots import TimeSeriesPlot
from utils import get_cable_informations
from fetch_data import clean_B_fields

from cable import SCUBASModel


folder="simulation/May2024/"
base_folder="dataset/May2024/"

def fetch_MAGE_simulations(stations=["frd", "stj", "had"]):
    global folder, base_folder
    for stn in stations:
        file = f"{base_folder}{stn.upper()}_May2024_MAGE.csv"
        d = pd.read_csv(file)
        d.rename(columns=dict(Var1="X",Var2="Y",Var3="Z"), inplace=True)
        d["Date"] = [dt.datetime(2024,5,9,20) + dt.timedelta(minutes=i) for i in range(len(d))]
        d = d[d.Date<dt.datetime(2024,5,12)]
        d.to_csv(f"{folder}{stn.lower()}_May2024_MAGE.csv", index=False, header=True, float_format="%g")
    return

def create_mage_stack_plot(dates=[dt.datetime(2024,5,10,12), dt.datetime(2024,5,12)]):
    global folder, base_folder
    import matplotlib.dates as mdates
    ts = TimeSeriesPlot(
        dates, f"", num_subplots=3, text_size=15,
        major_locator=mdates.HourLocator(byhour=range(0, 24, 4)),
        minor_locator=mdates.HourLocator(byhour=range(0, 24, 1)),
    )
    files = glob.glob(f"{folder}/*May2024_MAGE.csv")
    files.sort()
    od = {}
    print(files)
    for i, f in enumerate(files):
        d = pd.read_csv(f, parse_dates=["Date"])
        d.rename(columns=dict(Date="dates"), inplace=True)
        d.set_index("dates", inplace=True)
        d = d.interpolate(method="linear")
        d.reset_index(inplace=True)
        d = d.rename(columns=dict(dbn_nez="Y", dbe_nez="X", dbz_nez="Z"))
        od[f.split("/")[-1].split("_")[0]] = d
        logger.info(f"File {f}/{len(d)}")
        print(d.head())
    ts.add_mag(od, stations=["frd"])
    ts.add_mag(od, stations=["stj"])
    ts.add_mag(od, stations=["had"])
    ts.save("figures/MAGE.stack.png")
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

    fname, file = (
        f"{folder}{stn}_May2024_MAGE.csv",
        f"{base_folder}{stn}20240510psec.sec.csv"
    )
    df = pd.read_csv(fname, parse_dates=["Date"])
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
    ax = ts.add_curve(df.Date, df.Y,)
    ts.add_curve(o.Date, o.Y, ax = ax, color="b")
    ts.add_curve(data.tval, -data.N_nez, ax = ax, color="r", ylabel=r"$N_{nez}(Y)$ [nT]", ylim=[-1000, 1000])
    ax = ts.add_curve(df.Date, df.X,)
    ts.add_curve(o.Date, o.X, ax = ax, color="b")
    ts.add_curve(data.tval, data.E_nez, ax = ax, color="r", ylabel=r"$E_{nez}(X)$ [nT]", ylim=[-1000, 1000])
    ts.save(f"figures/MAGE.validation.{stn}.png")
    return

def run_May2024_storm():
    dates = [dt.datetime(2024,5,10,12), dt.datetime(2024,5,12)]
    segment_files = [
        ["simulation/May2024/frd_May2024_MAGE.csv"],
        ["simulation/May2024/frd_May2024_MAGE.csv"],
        ["simulation/May2024/stj_May2024_MAGE.csv"],
        ["simulation/May2024/stj_May2024_MAGE.csv"],
        ["simulation/May2024/stj_May2024_MAGE.csv"],
        ["simulation/May2024/stj_May2024_MAGE.csv"],
        ["simulation/May2024/stj_May2024_MAGE.csv"], 
        ["simulation/May2024/had_May2024_MAGE.csv"],
        ["simulation/May2024/had_May2024_MAGE.csv"],
    ]
    model = SCUBASModel(segment_files=segment_files)
    model.initialize_TL()
    model.run_cable_segment()
    data = model.cable.tot_params.reset_index().copy()
    data.to_csv("simulation/May2024/SCUBAS-Simulation-MAGE.csv", float_format="%g", index=False)
    return

def SuperMAG_compare_plots():
    smag_simulations = pd.read_csv("simulation/May2024/SCUBAS-Simulation-SuperMAG-Fitted.csv", parse_dates=["Time"]).set_index("Time")
    print(smag_simulations.head())
    scuba_simulations = pd.read_csv("simulation/May2024/SCUBAS-Simulation-3-Stations.csv", parse_dates=["Time"]).set_index("Time")
    mage_simulations = pd.read_csv("simulation/May2024/SCUBAS-Simulation-MAGE.csv", parse_dates=["Time"]).set_index("Time")
    print(mage_simulations.head())
    dates = [dt.datetime(2024,5,10,12), dt.datetime(2024,5,12)]
    import matplotlib.dates as mdates
    ts = TimeSeriesPlot(
        dates,
        major_locator=mdates.HourLocator(byhour=range(0, 24, 6)),
        minor_locator=mdates.HourLocator(byhour=range(0, 24, 3)),
        fig_title="Compare Voltages (R: 3 Stn, Bk: SMAG Fit, Bl: 3 Stn MAGE)", 
        text_size=15,
        num_subplots=1,
    )
    ax = ts.add_voltage(smag_simulations, xlabel="")
    ts.add_voltage(mage_simulations, ax=ax, color="b", xlabel="")
    ts.add_voltage(scuba_simulations, ax=ax, color="r", xlabel="Hours Since 10 May 2024, 12 UT")
    ts.save("simulation/May2024/SCUBAS-Simulation-Compare.png")
    ts.close()
    return

if __name__ == "__main__":
    # 1. Fetch o, dates based on MAGE simulations 10 and 11 th May (D)
    fetch_MAGE_simulations()
    # fetch_data_by_station("had")

    # TODO: Plot all the supermag dataset for different segments TS plot
    create_mage_stack_plot()
    stn_validation_plots("frd")
    stn_validation_plots("stj")
    stn_validation_plots("had")

    # TODO: Invoke scubas by the interpolated dataset 
    # run_May2024_storm()
    # SuperMAG_compare_plots()