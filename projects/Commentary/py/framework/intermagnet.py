import datetime as dt
import glob
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from loguru import logger


def get_tapering_function(t, p=0.1):
    """
    This method is resposible for generateing
    tapering function based on time sequence t
    and tapering coefficient p
    """
    T = len(t)
    P, P2 = int(T * p), int(T * p / 2)
    w = np.zeros_like(t)
    w[:P2] = 0.5 * (1 - np.cos(2 * np.pi * t[:P2] / P))
    w[P2 : T - P2] = 1.0
    w[T - P2 :] = 0.5 * (1 - np.cos(2 * np.pi * (t[-1] - t[T - P2 :]) / P))
    return w


@dataclass
class Bfield:
    """Class to represent the magnetic field vector."""

    bx: Sequence[float]
    by: Sequence[float]
    bz: Sequence[float]

    bh: Optional[Sequence[float]] = None
    bhang: Optional[Sequence[float]] = None
    date: Optional[Sequence[dt.datetime]] = None
    del_ta: Optional[float] = None
    dbx: Optional[Sequence[float]] = None
    dby: Optional[Sequence[float]] = None
    dbz: Optional[Sequence[float]] = None
    dbh: Optional[Sequence[float]] = None

    def __post_init__(self):
        self.bx = np.array(self.bx)
        self.by = np.array(self.by)
        self.bz = np.array(self.bz)
        self.bh = np.array(np.sqrt(self.bx**2 + self.by**2))
        self.bhang = np.array(np.degrees(np.pi + np.arctan2(self.by, self.bx)))

        self.dbx = np.array(np.gradient(self.bx, self.del_ta))
        self.dby = np.array(np.gradient(self.by, self.del_ta))
        self.dbz = np.array(np.gradient(self.bz, self.del_ta))
        self.dbh = np.array(np.gradient(self.bh, self.del_ta))
        self.dbhang = np.array(np.degrees(np.pi + np.arctan2(self.dby, self.dbx)))
        return

    def __len__(self):
        return len(self.bx)

    def __getitem__(self, item):
        return (self.bx[item], self.by[item], self.bz[item], self.bh[item])

    def __iter__(self):
        return iter((self.bx, self.by, self.bz, self.bh))

    def __repr__(self, item):
        return f"Bfield(bx={self.bx[item]}, by={self.by[item]}, bz={self.bz[item]}, bh={self.bh[item]})"

    @staticmethod
    def create_benchmark_bfield(
        folder: str = "data/benchmark/",
        files: Optional[Sequence[str]] = None,
    ) -> "Bfield":
        """Create a benchmark Bfield object."""
        files = glob.glob(f"{folder}*.csv")
        files.sort()
        logger.info(f"Found {len(files)} benchmark files in {folder}")
        if not files:
            logger.error(f"No benchmark files found in {folder}")
            raise FileNotFoundError("No benchmark files found in data/benchmark/")
        scale = (0.001 * np.exp(0.115 * 60)) / (0.001 * np.exp(0.115 * 55))
        DS = pd.concat([pd.read_csv(f, parse_dates=["datetime"]) for f in files])
        for key in ["x", "y", "z"]:
            DS[key + "_o"] = DS[key] * scale
        del_ta = (DS.datetime.iloc[1] - DS.datetime.iloc[0]).total_seconds()
        DS.x, DS.y, DS.z = (
            (DS.x - np.mean(DS.x.iloc[:60]))
            * scale
            * get_tapering_function(np.arange(len(DS)) * del_ta),
            (DS.y - np.mean(DS.y.iloc[:60]))
            * scale
            * get_tapering_function(np.arange(len(DS)) * del_ta),
            (DS.z - np.mean(DS.z.iloc[:60]))
            * scale
            * get_tapering_function(np.arange(len(DS)) * del_ta),
        )
        logger.info(f"Time step: {del_ta} seconds")
        bf = Bfield(
            bx=DS["x"].values,
            by=DS["y"].values,
            bz=DS["z"].values,
            date=DS["datetime"].dt.to_pydatetime(),
            del_ta=del_ta,
        )
        DS["Bh"] = np.sqrt(DS.x_o**2 + DS.y_o**2)
        DS["Bh"] = DS["Bh"] - np.mean(DS["Bh"].iloc[:60])
        setattr(bf, "bmag", bf.bh)
        setattr(bf, "bh", DS["Bh"].values)

        setattr(bf, "benchmark", True)
        setattr(bf, "benchmark_folder", folder)
        setattr(bf, "benchmark_files", files)
        setattr(bf, "benchmark_scale", scale)
        setattr(bf, "benchmark_dataset", DS)
        return bf

    @staticmethod
    def create_simulated_bfield(
        Am: Sequence = [200, 90, 30, 17, 8, 3.5],
        Tm_min: Sequence = [180, 80, 34, 15, 6, 3],
        Phim: Sequence = [10, 20, 30, 40, 50, 60],
        T_hours: float = 72,
        del_ta: float = 60.0,
        set_axis: str = "x",
    ) -> "Bfield":
        """Create a simulated Bfield object."""
        # Simulated data generation logic goes here
        # For now, we will just return a dummy Bfield object
        dta = int(T_hours * 3600 / del_ta)
        t = np.linspace(0, dta, dta, endpoint=False)
        b = np.zeros_like(t)
        for A, Phi, T in zip(Am, Phim, Tm_min):
            b += A * np.sin((2 * np.pi * t / T) + np.deg2rad(Phi))
        bf = Bfield(
            bx=b if set_axis == "x" else np.zeros_like(t),
            by=b if set_axis == "y" else np.zeros_like(t),
            bz=np.zeros_like(t),
            del_ta=del_ta,
        )
        return bf


@dataclass
class Efield:
    """Class to represent the electric field vector."""

    ex: Sequence[float]
    ey: Sequence[float]

    ez: Optional[Sequence[float]] = None
    eh: Optional[Sequence[float]] = None
    ehang: Optional[Sequence[float]] = None
    date: Optional[Sequence[dt.datetime]] = None
    del_ta: Optional[float] = None

    def __post_init__(self):
        self.ex = np.array(self.ex)
        self.ey = np.array(self.ey)
        self.ez = np.array(self.ez)
        self.eh = np.array(np.sqrt(self.ex**2 + self.ey**2))
        self.ehang = np.array(np.degrees(np.pi + np.arctan2(self.ey, self.ex)))
        return

    def __len__(self):
        return len(self.ex)

    def __getitem__(self, item):
        return (self.ex[item], self.ey[item], self.ez[item], self.eh[item])

    def __iter__(self):
        return iter((self.ex, self.ey, self.ez, self.eh))

    def __repr__(self, item):
        return f"Efield(ex={self.ex[item]}, ey={self.ey[item]}, ez={self.ez[item]}, eh={self.eh[item]})"
