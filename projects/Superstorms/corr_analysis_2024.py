import pandas as pd
import os
import datetime as dt
from loguru import logger
import numpy as np
from scubas.datasets import PROFILES

from mit import analog

from event_runs import (
    read_bfield_frames, 
    compile_cable_to_calculate_parameters,
    initialize_mag_fold,

)
from event_plots import (
    plot_correlation_function,
    create_correlation_function
)

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
        & (omni.time<dates[1])
    ]
    return omni



if __name__ == "__main__":
    initialize_mag_fold()
    FRD, HAD, STJ = (
        [".scubas_config/frd20240510psec.sec.txt"],
        [".scubas_config/had20240510psec.sec.txt"], 
        [".scubas_config/stj20240510vsec.sec.txt"]
    )
    stns, dates = (
        ["FRD", "STJ", "HAD"], 
        [dt.datetime(2024,5,10,12), dt.datetime(2024,5,12)]
    )
    frames = read_bfield_frames(stns, [FRD, HAD, STJ])

    profiles = [
        PROFILES.CS_W, PROFILES.DO_1, PROFILES.DO_2, PROFILES.DO_3,
        PROFILES.DO_4, PROFILES.DO_5, PROFILES.MAR, PROFILES.DO_6,
        PROFILES.CS_E
    ]
    tlines, cable = compile_cable_to_calculate_parameters(FRD, STJ, HAD, profiles)
    cable.tot_params = cable.tot_params[
        (cable.tot_params.index>=dates[0])
        & (cable.tot_params.index<=dates[1])
    ]
    omni = _load_omni_(dates)
    for c in cable.tot_params.columns:
        omni[c] = np.array(cable.tot_params[c].tolist()[::60])
    
    rolling_corr = dict(
        time=omni.time
    )
    print(omni.columns)
    X_params, Y_params, ylabels = (
        [
            "B", "ProtonDensity", "BzGSM", "AE"
        ], [
            "V(v)",
        ], [
            r"B, nT", r"n, /cc", r"$B_z$", r"$A_E$"
        ]
    )
    # window = 10
    # for x in X_params:
    #     for y in Y_params:
    #         rolling_corr[f"{x}_{y}"] = (
    #             omni[x].interpolate(method='quadratic').rolling(window=window).corr(omni[y])
    #         )
    
    # plot_correlation_function(
    #     rolling_corr, 
    #     X_params, Y_params, ylabels=ylabels,
    #     xlim=[dt.datetime(2024,5,10,16), dt.datetime(2024,5,10,18)], 
    # )
    omni = omni[
        (omni.time>=dt.datetime(2024,5,10,17))
        & (omni.time<=dt.datetime(2024,5,10,18))
    ] 
    omni = omni.interpolate()
    for x in X_params:
        print(omni[x].isna().any())
        print(analog.joint_differential_entropy(
            omni["V(v)"].interpolate(method='quadratic'),
            omni[x].interpolate(method='quadratic'),
        ))
    create_correlation_function(
        omni[
            (omni.time>=dt.datetime(2024,5,10,17))
            & (omni.time<=dt.datetime(2024,5,10,18))
        ], X_params, ylabels
    )