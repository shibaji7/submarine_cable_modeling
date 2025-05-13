from types import SimpleNamespace

import numpy as np
import pandas as pd
from fetch_data import read_iaga
from loguru import logger
from scubas.utils import GreatCircle


def parse_csv(fname="datasets/20250211-20-39-supermag.csv"):
    o = pd.read_csv(fname, parse_dates=["Date_UTC"])
    iagas = o.IAGA.unique()
    for iaga in iagas:
        logger.info(f"IAGA: {iaga}")
        x = o[o.IAGA == iaga][
            [
                "Date_UTC",
                "dbn_geo",
                "dbe_geo",
                "dbz_geo",
                "IAGA",
                "GEOLON",
                "GEOLAT",
            ]
        ]
        x["F"] = np.sqrt(x.dbn_geo**2 + x.dbe_geo**2 + x.dbz_geo**2)
        x.rename(
            columns={
                "Date_UTC": "Date",
                "dbn_geo": "Y",
                "dbe_geo": "X",
                "dbz_geo": "Z",
                "IAGA": "code",
                "GEOLON": "lon",
                "GEOLAT": "lat",
            },
            inplace=True,
        )
        x.to_csv(f"datasets/{iaga}.csv", header=True, index=False)
    return


def get_nearest_station(lat, lon, codes=["CNB", "CTA", "GUA", "KAK"]):
    distances = []
    for code in codes:
        # fname = f"datasets/{code.upper()}.csv"
        fname = f"datasets/IM.{code.lower()}.min.txt"
        o, header = read_iaga(fname, return_header=True)
        glat, glon = (
            float(header["geodetic latitude"]),
            float(header["geodetic longitude"]),
        )
        gc = GreatCircle(
            SimpleNamespace(**dict(lat=lat, lon=lon)),
            SimpleNamespace(**dict(lat=glat, lon=glon)),
        )
        distances.append(gc.great_circle())
    code = codes[np.argmin(distances)]
    return code, f"datasets/IM.{code.lower()}.min.csv"


def load_modify_file(code):
    fname = f"datasets/{code.upper()}.csv"
    o = pd.read_csv(fname, parse_dates=["Date"])
    return


if __name__ == "__main__":
    parse_csv()
