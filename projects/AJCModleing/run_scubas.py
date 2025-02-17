import sys
import datetime as dt
# sys.path.append("py/")
from fetch_data import (
    clean_B_fields, load_speadas, _load_omni_
)
from cable import SCUBASModel
from fetch_cable_routes import (
    get_cable_route,
    calculate_bathymetry_byLITHO1,
    calculate_conductive_profiles
)

import os

def run_May2024_storm():
    os.makedirs("simulation/May2024/", exist_ok=True)
    dates = [dt.datetime(2024,5,10,12), dt.datetime(2024,5,12)]
    themis_fgm, themis_mom = load_speadas(dates)

    # Create dSegments and Profiles
    # 1. From cable routes run following methods
    #   o = get_cable_route()
    #   d = calculate_bathymetry_byLITHO1(o)
    #   profiles = calculate_conductive_profiles(d)
    # This will help to check how d, profiles will help computing segment files
    o = get_cable_route()
    d = calculate_bathymetry_byLITHO1(o)
    print(d) 
    # 2. Manually check and assign the segment files to each profile
    # Probably remove some points in 'd' and then send for profilings
    profiles = calculate_conductive_profiles(d)
    
    # 3. Conduuct a mag data cleaining by invoking  clean_B_fields 
    #       (will make sure files are readable)

    # 4. Run SCUBASModel for specific AJC cable
    #   Note: Make sure you have modified get_cable_informations 

    # 5. Run the folloding lines to finally run the code and save (plot).

    # segment_files = [
    #     ["dataset/May2024/frd20240510psec.sec.csv"],
    #     ["dataset/May2024/frd20240510psec.sec.csv"],
    #     ["dataset/May2024/stj20240510psec.sec.csv"],
    #     ["dataset/May2024/stj20240510psec.sec.csv"],
    #     ["dataset/May2024/stj20240510psec.sec.csv"],
    #     ["dataset/May2024/stj20240510psec.sec.csv"],
    #     ["dataset/May2024/stj20240510psec.sec.csv"], 
    #     ["dataset/May2024/had20240510psec.sec.csv"],
    #     ["dataset/May2024/had20240510psec.sec.csv"],
    # ]
    # clean_B_fields(
    #     [
    #         "FRD", "STJ", "HAD"
    #     ], [
    #         ["dataset/May2024/frd20240510psec.sec.txt"],
    #         ["dataset/May2024/stj20240510psec.sec.txt"],
    #         ["dataset/May2024/had20240510psec.sec.txt"]
    #     ]
    # )
    # model = SCUBASModel(segment_files=segment_files)
    # model.initialize_TL()
    # model.run_cable_segment()
    # data = model.cable.tot_params.reset_index().copy()
    # data.to_csv("simulation/May2024/SCUBAS-Simulation-3-Stations.csv", float_format="%g", index=False)
    return

if __name__ == "__main__":
    run_May2024_storm()