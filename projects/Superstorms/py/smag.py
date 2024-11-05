import xarray as xr
import datetime as dt
from loguru import logger
import numpy as np
import aacgmv2
import pandas as pd
from types import SimpleNamespace
import os

from event_plots import TimeSeriesPlot

def get_cable_informations(kind="TAT-8"):
    if kind == "TAT-8":
        cable = SimpleNamespace(**dict(
            cable_seg = [
                dict(
                    initial=dict(lat=39.6, lon=-74.33), 
                    final=dict(lat=38.79, lon=-72.62)
                ),
                dict(
                    initial=dict(lat=38.79, lon=-72.62), 
                    final=dict(lat=37.11, lon=-68.94)
                ),
                dict(
                    initial=dict(lat=37.11, lon=-68.94), 
                    final=dict(lat=39.80, lon=-48.20)
                ),
                dict(
                    initial=dict(lat=39.80, lon=-48.20), 
                    final=dict(lat=40.81, lon=-45.19)
                ),
                dict(
                    initial=dict(lat=40.81, lon=-45.19), 
                    final=dict(lat=43.15, lon=-39.16)
                ),
                dict(
                    initial=dict(lat=43.15, lon=-39.16), 
                    final=dict(lat=44.83, lon=-34.48)
                ),
                dict(
                    initial=dict(lat=44.83, lon=-34.48), 
                    final=dict(lat=46.51, lon=-22.43)
                ),
                dict(
                    initial=dict(lat=46.51, lon=-22.43), 
                    final=dict(lat=47.85, lon=-9.05)
                ),
                dict(
                    initial=dict(lat=47.85, lon=-9.05), 
                    final=dict(lat=50.79, lon=-4.55)
                ),
            ]
        ))
    for seg in cable.cable_seg:
        seg["center"] = dict(
            lat=0.5*(seg["initial"]["lat"]+seg["final"]["lat"]),
            lon=0.5*(seg["initial"]["lon"]+seg["final"]["lon"]),
        )
    return cable

def read_dataset(date):
    file = f".scubas_config/{date.strftime('%Y%m%d')}.north.schavec-mlt-supermag.60s.rev-0006.ncdf"
    obj = xr.open_dataset(file)
    ddates, dates = [], []
    for y,m,d,H,M,S in zip(
            obj.variables["time_yr"].values.astype(int),
            obj.variables["time_mo"].values.astype(int),
            obj.variables["time_dy"].values.astype(int),
            obj.variables["time_hr"].values.astype(int),
            obj.variables["time_mt"].values.astype(int),
            obj.variables["time_sc"].values.astype(int),
        ):
        ddates.extend([dt.datetime(y,m,d,H,M,S)]*600)
        dates.append(dt.datetime(y,m,d,H,M,S))
    o = pd.DataFrame()
    o["dates"] = ddates
    o["mlat"] = obj.variables["mlat"].values.ravel()
    o["mlt"] = obj.variables["mlt"].values.ravel()
    o["mlon"] = np.mod(obj.variables["mlon"].values.ravel()+180, 360) - 180
    o["dbn_nez"] = obj.variables["dbn_nez"].values.ravel()
    o["dbe_nez"] = obj.variables["dbe_nez"].values.ravel()
    o["dbz_nez"] = obj.variables["dbz_nez"].values.ravel()
    return o, dates

def get_dataset_from_netcdf(date, glat, glon, df, dmlat=2, dmlt=1):
    mlat, mlon, mlt = aacgmv2.get_aacgm_coord(glat, glon, 300, date)
    logger.info(f"Converted mlat, mlon, mlt: {mlat}, {mlon}, {mlt}")
    x = df[
        (df.dates == date)
        & (df.mlt>=np.mod(mlt-dmlt, 24))
        & (df.mlt<=np.mod(mlt+dmlt, 24))
        & (df.mlat>=mlat-dmlat)
        & (df.mlat<=mlat+dmlat)
    ]
    logger.info(f"Number of data points againts this bin {len(x)}")
    x = x.groupby("dates").mean().reset_index()
    return x

def fetch_dataset_from_netcdf_by_date_range(dates, cable, df, dmlat=2, dmlt=1):
    for i, seg in enumerate(cable.cable_seg):
        glat, glon = seg["center"]["lat"], seg["center"]["lon"]
        frame = pd.DataFrame()
        for date in dates:
            x = get_dataset_from_netcdf(date, glat, glon, df, dmlat, dmlt)
            frame = pd.concat([frame, x])
        
        fname = f"simulation/May2024/cable_segment_{i}.csv"
        frame.to_csv(fname, float_format="%g", index=False, header=True)
    return

def interpolate_smag():
    os.makedirs("simulation/May2024/", exist_ok=True)
    cable = get_cable_informations()
    O, Dates = pd.DataFrame(), []
    for d in [dt.datetime(2024,5,10), dt.datetime(2024,5,11)]:
        o, dates = read_dataset(d)
        O = pd.concat([O, o])
        Dates.extend(dates)
    fetch_dataset_from_netcdf_by_date_range(
        Dates, cable, O, 
        dmlat=2, dmlt=1
    )
    return

if __name__ == "__main__":
    # 1. Fetch o, dates based on 10 and 11 th May (D)
    # 2. Call fetch_dataset_from_netcdf_by_date_range for whole date range 
    #       and save to local files. (D)
    # TODO: 3. May be plan to add 12th May
    interpolate_smag()

    # TODO: Plot all the supermag dataset for 9 different segments TS plot
    
    # TODO: Invoke scubas by the interpolated dataset 