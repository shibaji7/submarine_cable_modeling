import datetime as dt
import glob
import sys

import matplotlib.dates as mdates
import pandas as pd

sys.path.append("py/")
from fetch_data import _load_omni_
from plots import TimeSeriesPlot

import numpy as np
from scipy.signal import butter, filtfilt
def filter_data(data, cutoff, fs, order=4):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    
    # Get the filter coefficients (b, a)
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    
    # Apply the filter with zero phase offset
    y = filtfilt(b, a, data)
    return y


def read_files_as_pandas(folder):
    df = []
    print(f"Folder: {folder}")
    files = glob.glob(folder + "/*.csv")
    files.sort()
    for f in files:
        print(f"\tFile: {f}")
        o = pd.read_csv(f, parse_dates=["Timestamp"])
        df.append(o)
    df = pd.concat(df)
    return df


def plot_datasets_instack(
    dfs,
    dates=[dt.datetime(2024, 5, 10, 12), dt.datetime(2024, 5, 12)],
    fig_title="",
    volt_key="V_(S1)_PAD (1S) (V)",
    curr_key="I_(S1)_PAD (1S) (mA)",
    fname="figures/AJC.png",
    is_filter=False,
):
    omni = _load_omni_(dates)
    ts = TimeSeriesPlot(
        dates,
        fig_title,
        num_subplots=len(dfs),
        text_size=15,
        major_locator=mdates.HourLocator(byhour=range(0, 24, 4)),
    )
    for _, df in enumerate(dfs):
        df.Timestamp = df.Timestamp + dt.timedelta(minutes=15)
        # df = df[
        #     (df.Timestamp>=dates[0])
        #     & (df.Timestamp<=dates[-1])
        # ]
        ax = ts._add_axis()
        ax.plot(df.Timestamp, df[volt_key], "b.", ls="None", ms=0.5, alpha=0.7)
        ax.set_xlabel("Time, UT")
        ax.set_ylabel("Voltage, V", fontdict=dict(color="b"))
        tax = ax.twinx()
        tax.set_ylabel("IMF [B], nT", fontdict=dict(color="r"))
        tax.plot(omni.time, omni.B, "r.", ls="None", ms=0.5, alpha=1.0)
        if is_filter:
            df = df.set_index('Timestamp')
            df['interpolated'] = df[volt_key].interpolate(method='time')
            data = np.array(df["interpolated"].tolist())
            fs = 1/(df.index.tolist()[1]-df.index.tolist()[0]).total_seconds()
            print(fs, df.index.tolist()[1], df.index.tolist()[0], 
                (df.index.tolist()[1]-df.index.tolist()[0]).total_seconds())
            filtered_data = filter_data(data, cutoff=1/(0.5*3600), fs=fs, order=6)
            ax.plot(df.index, filtered_data, "k", ls="-", ms=0.5, alpha=0.7)
    ts.save(fname)

    print(omni.columns)
    return


if __name__ == "__main__":
    # fold_name = "Oxford Falls (AUSTRALIA)"
    # df = read_files_as_pandas(f"dataset/Operators/AJC PFE DATA/{fold_name}")
    # stn, station_code, data_cad = "OXF", 2, 1
    # plot_datasets_instack(
    #     [df],
    #     fig_title=f"Date: 10-12 May 2024; Stn: {stn.lower()}/aus",
    #     volt_key=f"V_(S{station_code})_{stn} ({data_cad}S) (V)",
    #     curr_key=f"I_(S{station_code})_{stn} ({data_cad}S) (mA)",
    #     fname=f"figures/ajc_{stn.lower()}.png",
    # )

    # fold_name = "Tumon Bay (GUAM)"
    # df = read_files_as_pandas(f"dataset/Operators/AJC PFE DATA/{fold_name}")
    # stn, station_code, data_cad = "TMB", 4, 1
    # plot_datasets_instack(
    #     [df],
    #     fig_title=f"Date: 10-12 May 2024; Stn: {stn.lower()}/aus",
    #     volt_key=f"S_4_V (V)",
    #     curr_key=f"S_4_I (mA)",
    #     fname=f"figures/ajc_{stn.lower()}.png",
    # )

    fold_name = "Maruyama (JAPAN)"
    df = read_files_as_pandas(f"dataset/Operators/AJC PFE DATA/{fold_name}")
    stn, station_code, data_cad = "mar", 4, 1
    plot_datasets_instack(
        [df],
        fig_title=f"Date: 10-12 May 2024; Stn: {stn.lower()}/aus",
        volt_key=f"V_(S10)_MRU (5S) (V)",
        curr_key=f"I_(S10)_MRU(5S) (mA)",
        fname=f"figures/ajc_{stn.lower()}.png",
        is_filter=True,
    )

    # fold_name = "Tanguisson (Guam)"
    fold_name = "Tanguisson (GUAM)"
    df = read_files_as_pandas(f"dataset/Operators/AJC PFE DATA/{fold_name}")
    stn, station_code, data_cad = "tng", 4, 1
    plot_datasets_instack(
        [df],
        fig_title=f"Date: 10-12 May 2024; Stn: {stn.lower()}/guam",
        volt_key=f"V_S6_TNG (V)",
        curr_key=f"I_S6_TNG (mA)",
        fname=f"figures/ajc_{stn.lower()}.png",
        is_filter=True,
    )
