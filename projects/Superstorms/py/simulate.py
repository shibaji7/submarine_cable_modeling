import datetime as dt
import sys

sys.path.append("py/")
import os

from cable import SCUBASModel
from fetch_data import clean_B_fields, load_speadas


def run_May2024_storm():
    os.makedirs("simulation/May2024/", exist_ok=True)
    dates = [dt.datetime(2024, 5, 10, 12), dt.datetime(2024, 5, 12)]
    themis_fgm, themis_mom = load_speadas(dates)
    segment_files = [
        ["dataset/May2024/frd20240510psec.sec.csv"],
        ["dataset/May2024/frd20240510psec.sec.csv"],
        ["dataset/May2024/stj20240510psec.sec.csv"],
        ["dataset/May2024/stj20240510psec.sec.csv"],
        ["dataset/May2024/stj20240510psec.sec.csv"],
        ["dataset/May2024/stj20240510psec.sec.csv"],
        ["dataset/May2024/stj20240510psec.sec.csv"],
        ["dataset/May2024/had20240510psec.sec.csv"],
        ["dataset/May2024/had20240510psec.sec.csv"],
    ]
    clean_B_fields(
        ["FRD", "STJ", "HAD"],
        [
            ["dataset/May2024/frd20240510psec.sec.txt"],
            ["dataset/May2024/stj20240510psec.sec.txt"],
            ["dataset/May2024/had20240510psec.sec.txt"],
        ],
    )
    model = SCUBASModel(segment_files=segment_files)
    model.initialize_TL()
    model.run_cable_segment()
    data = model.cable.tot_params.reset_index().copy()
    data.to_csv(
        "simulation/May2024/SCUBAS-Simulation-3-Stations.csv",
        float_format="%g",
        index=False,
    )
    return


if __name__ == "__main__":
    run_May2024_storm()
