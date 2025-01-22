import pandas as pd
import glob
import datetime as dt
import sys
sys.path.append("py/")
from plots import TimeSeriesPlot

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
        dfs, dates=[dt.datetime(2024,5,8), dt.datetime(2024,5,13)], 
        fig_title="", volt_key="V_(S1)_PAD (1S) (V)", curr_key="I_(S1)_PAD (1S) (mA)",
        fname="figures/AJC.png"
    ):
    ts = TimeSeriesPlot(
        dates, fig_title, num_subplots=len(dfs), text_size=15
    )
    for _, df in enumerate(dfs):
        ax = ts._add_axis()
        print(df.head())
        ax.plot(df.Timestamp, df[volt_key], "b.", ls="None", ms=0.5, alpha=0.7)
        ax.set_xlabel("Time, UT")
        ax.set_ylabel("Voltage, V", fontdict=dict(color="b"))
        ax = ax.twinx()
        ax.set_ylabel("Current, mA", fontdict=dict(color="r"))
        ax.plot(df.Timestamp, df[curr_key], "r.", ls="None", ms=0.5, alpha=0.7)
    ts.save(fname)
    return


if __name__ == "__main__":
    df = read_files_as_pandas(
        "dataset/Operators/AJC PFE DATA/Shima (JAPAN)"
    )
    plot_datasets_instack(
        [df], fig_title="Date: 09-13 May 2024; Stn: shi/jpn",
        volt_key="V_(S10)_SHI (5S) (V)", curr_key="I_(S10)_SHI (5S) (mA)",
        fname="figures/ajc_shi.png"
    )