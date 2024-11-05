import sys
sys.path.append("py/")
from fetch_data import clean_B_fields
from cable import SCUBASModel

if __name__ == "__main__":
    segment_files = [
        ["dataset/May2024/frd20240510psec.sec.csv"],
        ["dataset/May2024/frd20240510psec.sec.csv"],
        ["dataset/May2024/stj20240510vsec.sec.csv"],
        ["dataset/May2024/stj20240510vsec.sec.csv"],
        ["dataset/May2024/stj20240510vsec.sec.csv"],
        ["dataset/May2024/stj20240510vsec.sec.csv"],
        ["dataset/May2024/stj20240510vsec.sec.csv"], 
        ["dataset/May2024/had20240510psec.sec.csv"],
        ["dataset/May2024/had20240510psec.sec.csv"],
    ]
    clean_B_fields(
        [
            "FRD", "STJ", "HAD"
        ], [
            ["dataset/May2024/frd20240510psec.sec.txt"],
            ["dataset/May2024/stj20240510vsec.sec.txt"],
            ["dataset/May2024/had20240510psec.sec.txt"]
        ]
    )
    model = SCUBASModel(segment_files=segment_files)
    model.initialize_TL()