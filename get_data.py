"""get_data.py: Module is used to implement download magnetometer data, store, read and do analysis FFT B(f)"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib.colors import LogNorm

import os
os.system("rm -rf `find -type d -name .ipynb_checkpoints`:")
import numpy as np
import datetime as dt
import pandas as pd
from supermag_api import *

UID = "shibaji7"
def download_magnetometer_data_by_stations_dates(start, end, stns, local_dir = "tmp/", flags = "mlt,mag,geo,decl,sza", 
                                                coord="geo"):
    """
    Download and store the data in CSV format
    """
    if not os.path.exists(local_dir): os.system("mkdir " + local_dir)
    dur = int((end-start).total_seconds())
    frames = {}
    for stn in stns:
        fname = local_dir + "%s_%s_%s.csv"%(stn,start.strftime("%Y%m%dT%H%M"), end.strftime("%Y%m%dT%H%M"))
        if not os.path.exists(fname):
            status, sm_data = SuperMAGGetData(UID, start.strftime("%Y-%m-%dT%H%M"), dur, flags, stn, FORMAT="pandas")
            if status == 1:
                sm_data["time"] = sm_data.tval.apply(lambda x: dt.datetime.utcfromtimestamp(x))
                sm_data.N = sm_data.N.apply(lambda x: x[coord])
                sm_data.E = sm_data.E.apply(lambda x: x[coord])
                sm_data.Z = sm_data.Z.apply(lambda x: x[coord])
                sm_data.to_csv(fname, index=False, header=True)
                frames[stn] = sm_data
        else: frames[stn] = pd.read_csv(fname, parse_dates=["time"])
    return frames

def plotB(frames, kind="nez", ylim=[-100,100]):
    """
    Plot B(t): Three components of B field
    """
    fig, axes = plt.subplots(figsize=(6, 3*len(frames.keys())), dpi=150, nrows=len(frames.keys()), 
                             ncols=1, sharex="all", sharey="all")
    if len(frames.keys()) == 1: axes = [axes]
    for k, ax in zip(frames.keys(), axes):
        frame = frames[k]
        ax.xaxis.set_major_formatter(DateFormatter(r"$%H^{%M}$"))
        ax.plot(frame.time, frame.N, "r", ls="-", lw=0.8, label=r"$\vec{B}_N$")
        ax.plot(frame.time, frame.E, "b", ls="-", lw=0.8, label=r"$\vec{B}_E$")
        ax.plot(frame.time, frame.Z, "k", ls="-", lw=0.8, label=r"$\vec{B}_Z$")
        ax.legend(loc=1)
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_xlim(frame.time.tolist()[0], frame.time.tolist()[-1])
        ax.set_ylabel(r"Magnetic Flux ($\vec{B}), nT$")
    ax.set_xlabel("Time, UT")
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    return

def Btime2Bfreq(frames, dw=30, sample_rate=1, component="N"):
    """
    FFT to convert B(t) -> B(f)
    frames: Dictionary of pandas
    dw: Window in minutes
    sample_rate: Sample rate in miutes
    """
    from scipy.fft import fft, fftfreq, fftshift
    Bf_frames = {}
    win, sr = dw*60, sample_rate
    N = win * sr
    for k in frames.keys():
        frame = frames[k]
        frame = frame.set_index("time").resample("1s").interpolate().reset_index()
        Bf_frames[k] = {"f":[], "t":[], "y":[], "component": component}
        time, Bt = frame.time.tolist(), frame[component].tolist()
        for ix in range(int(win/2), len(time)-int(win/2)):
            yf = fft(Bt[ix-int(win/2):ix+int(win/2)])
            xf = fftshift(fftfreq(N, 1 / sr))
            Bf_frames[k]["y"].append(yf[int(win/2):])
            Bf_frames[k]["f"].append(xf[int(win/2):])
            Bf_frames[k]["t"].append(int(win/2)*[time[ix]])
    return Bf_frames

def plotBf(frames, ylim=[1e-3, 1e0]):
    """
    Plot B(f): One components of B field
    """
    fig, axes = plt.subplots(figsize=(6, 3*len(frames.keys())), dpi=150, nrows=len(frames.keys()), 
                             ncols=1, sharex="all", sharey="all")
    if len(frames.keys()) == 1: axes = [axes]
    norm=LogNorm()
    for k, ax in zip(frames.keys(), axes):
        frame = frames[k]
        ax.xaxis.set_major_formatter(DateFormatter(r"$%H^{%M}$"))
        ax.xaxis_date()
        f, t, y = frame["f"], [ut[0] for ut in frame["t"]], np.array(frame["y"]).T
        ax.contourf(t, f[0], np.absolute(y), cmap=plt.cm.Spectral, norm=norm)
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_xlim(t[0], t[-1])
        ax.set_xlabel("Time, UT")
        ax.set_ylabel(r"$f_0$, Hz")
        ax.set_yscale("log")
        pos = ax.get_position()
        cpos = [pos.x1 + 0.025, pos.y0 + 0.0125,
                0.015, pos.height * 0.8]                # this list defines (left, bottom, width, height
        cax = fig.add_axes(cpos)
        cb2 = mpl.colorbar.ColorbarBase(cax, cmap=plt.cm.Spectral, norm=norm, spacing="uniform", orientation="vertical")
        cb2.set_label(r"PSD [$|B(f)|$], $nT.Hz^{-1}$")
        ax.text(0.01, 1.05, "Date: "+t[0].strftime("%Y-%m-%d"), ha="left", va="center", transform=ax.transAxes)
        ax.text(0.99, 1.05, "Component: "+r"$\vec{B}_{%s}$"%frame["component"], ha="right", va="center", transform=ax.transAxes)
    ax.set_xlabel("Time, UT")
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    return

if __name__ == "__main__":
    frames = download_magnetometer_data_by_stations_dates(dt.datetime(2017,9,6,12), dt.datetime(2017,9,6,13), ["STJ"])
    frames = download_magnetometer_data_by_stations_dates(dt.datetime(2017,9,6,11), dt.datetime(2017,9,6,14,1), ["STJ"])
    plotB(frames)
    Bf_frames = Btime2Bfreq(frames)
    plotBf(Bf_frames)
    pass