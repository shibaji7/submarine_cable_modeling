import sys
import datetime as dt
sys.path.append("py/")
from fetch_data import (
    clean_B_fields, load_speadas, _load_omni_
)
from cable import SCUBASModel

def run_May2024_storm():
    dates = [dt.datetime(2024,5,10,12), dt.datetime(2024,5,12)]
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
        [
            "FRD", "STJ", "HAD"
        ], [
            ["dataset/May2024/frd20240510psec.sec.txt"],
            ["dataset/May2024/stj20240510psec.sec.txt"],
            ["dataset/May2024/had20240510psec.sec.txt"]
        ]
    )
    model = SCUBASModel(segment_files=segment_files)
    model.initialize_TL()
    model.run_cable_segment()
    return

if __name__ == "__main__":
    run_May2024_storm()