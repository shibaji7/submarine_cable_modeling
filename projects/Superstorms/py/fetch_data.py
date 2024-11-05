import os
import pyspedas
from loguru import logger
import pandas as pd
import datetime as dt
os.environ["OMNIDATA_PATH"] = "/home/shibaji/omni/"

def _load_omni_(dates, res=1):
    import pyomnidata
    logger.info(f"OMNIDATA_PATH: {os.environ['OMNIDATA_PATH']}")
    pyomnidata.UpdateLocalData()
    omni = pd.DataFrame(
        pyomnidata.GetOMNI(dates[0].year,Res=res)
    )
    omni["time"] = omni.apply(
        lambda r: (
            dt.datetime(
                int(str(r.Date)[:4]), 
                int(str(r.Date)[4:6]),
                int(str(r.Date)[6:].replace(".0","")) 
            ) 
            + dt.timedelta(hours=r.ut)
        ), 
        axis=1
    )
    omni = omni[
        (omni.time>=dates[0])
        & (omni.time<=dates[1])
    ]
    return omni

def load_speadas(dates, probe="c"):
    time_range = [
        dates[0].strftime("%Y-%m-%d/%H:%M"),
        dates[1].strftime("%Y-%m-%d/%H:%M")
    ]
    data_fgm = pyspedas.themis.fgm(
        probe=probe, trange=time_range, 
        time_clip=True, no_update=False,notplot=True
    )
    data_mom = pyspedas.themis.mom(
        probe=probe, trange=time_range,
        notplot=True,no_update=False,time_clip=True
    )
    pdyn = {
        "x": data_mom["thc_peem_density"]["x"], 
        "y": data_mom["thc_peem_density"]["y"]*1.67*(10**(-6))*0.5*np.nansum(
            data_mom["thc_peim_velocity_gse"]["y"]**2, axis=1
        )
    }
    data_mom["pdyn"] = pdyn
    return data_fgm, data_mom