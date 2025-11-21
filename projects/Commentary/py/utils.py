import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import constants as C
from scubas.datasets import PROFILES, Site

sites = [
    Site.init(
        conductivities=[1 / 3, 1 / 3000],
        thicknesses=[10000, np.inf],
        names=[
            "Sediments",
            "Lower Mantle",
        ],
        description="This model is modified for David's cometary paper",
        site_name="Two layer Earth model",
    ),
    Site.init(
        conductivities=[1 / 3000, 1 / 3],
        thicknesses=[100000, np.inf],
        names=[
            "Sediments",
            "Lower Mantle",
        ],
        description="This model is modified for David's cometary paper",
        site_name="Two layer Earth model",
    ),
    Site.init(
        conductivities=[1 / 3],
        thicknesses=[np.inf],
        names=[
            "Lower Mantle",
        ],
        description="This model is modified for David's cometary paper",
        site_name="Unform Earth",
    ),
    Site.init(
        conductivities=[1 / 3000],
        thicknesses=[np.inf],
        names=[
            "Lower Mantle",
        ],
        description="This model is modified for David's cometary paper",
        site_name="Unform Earth",
    ),
    PROFILES.DB,
]


def run_benchmark_sim(snum=0, dates=[]):
    from scubas.utils import fft, ifft

    site = sites[snum]
    ds, del_t = get_benchmark_datasets()
    Bxf, f = fft(ds.x, del_t)
    Ets = dict(
        Y=-ifft(np.array(calcTFx(site, freqs=f).E2B) * Bxf),
    )
    ds["Ey"] = -ifft(np.array(calcTFx(site, freqs=f).E2B) * Bxf)
    Byf, f = fft(ds.y, del_t)
    Ets["X"] = ifft(np.array(calcTFx(site, freqs=f).E2B) * Byf)
    ds["Ex"] = ifft(np.array(calcTFx(site, freqs=f).E2B) * Byf)
    if len(dates) >= 2:
        ds = ds[(ds.datetime >= dates[0]) & (ds.datetime <= dates[1])]
    return ds, Ets, del_t


def get_nerc_dataset():
    datasets = pd.read_csv("datasets/nerc.csv")
    return datasets


def get_benchmark_datasets(a_scale=None):
    a_55, a_60 = 0.001 * np.exp(0.115 * 55), 0.001 * np.exp(0.115 * 60)
    a_scale = a_scale if a_scale else a_60 / a_55
    files = glob.glob("datasets/OTT*.csv")
    files.sort()
    datasets = pd.concat([pd.read_csv(f, parse_dates=["datetime"]) for f in files])
    for key in ["x", "y", "z"]:
        datasets[key + "_o"] = datasets[key] * a_scale
    datasets.x = datasets.x - np.mean(datasets.x.iloc[:60])
    datasets.y = datasets.y - np.mean(datasets.y.iloc[:60])
    datasets.z = datasets.z - np.mean(datasets.z.iloc[:60])
    del_t = (datasets.datetime.iloc[1] - datasets.datetime.iloc[0]).total_seconds()
    for key in ["x", "y", "z"]:
        datasets[key] = datasets[key] * get_tapering_function(
            np.arange(len(datasets)) * del_t
        )
    datasets["dx"] = np.diff(datasets.x, prepend=datasets.x.iloc[0]) / del_t
    datasets["dy"] = np.diff(datasets.y, prepend=datasets.y.iloc[0]) / del_t
    datasets["dz"] = np.diff(datasets.z, prepend=datasets.z.iloc[0]) / del_t
    return datasets, del_t


def calcZ(site, freqs):
    omega = 2 * C.pi * freqs
    sigma_s = 1 / site.layers[0].resistivity
    k2 = 1.0j * omega * C.mu_0 * sigma_s
    k = np.sqrt(k2)
    Zo = (1.0j * omega * C.mu_0 / k) / (C.mu_0 / 1.0e-3)
    return Zo


def calcTF(site, freqs):
    """
    Calculate the transfer functions.
    """
    Zo = calcZ(site, freqs)
    Zd = Zo
    omega = 2 * C.pi * freqs
    sigma_s = 1 / site.layers[0].resistivity
    k2 = 1.0j * omega * C.mu_0 * sigma_s
    k = np.sqrt(k2)
    kd = k * 0

    func = Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Zo))
    tf = pd.DataFrame()
    tf["freq"], tf["E2B"] = freqs, func
    return tf


def calcTFx(site, freqs):
    """
    Calculate the transfer functions.
    """
    Zo = calcZ(site, freqs)
    Zd = site.calcZ(freqs)[1]
    omega = 2 * C.pi * freqs
    sigma_s = 1 / site.layers[0].resistivity
    k2 = 1.0j * omega * C.mu_0 * sigma_s
    k = np.sqrt(k2)
    kd = k * 0

    func = Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Zo))
    tf = pd.DataFrame()
    tf["freq"], tf["E2B"] = freqs, func
    return tf


def abline():
    gca = plt.gca()
    gca.set_autoscale_on(False)
    gca.plot(gca.get_xlim(), gca.get_ylim(), ls="-.", lw=0.4, color="k")
    return


def abxline():
    gca = plt.gca()
    gca.set_autoscale_on(False)
    gca.plot(gca.get_xlim(), [4e-1, 4e-1], ls="-.", lw=0.4, color="k")
    return


def xfft(signal, sr=1):
    fft_signal = np.fft.fft(signal)
    frequencies = np.fft.fftfreq(len(signal), 1 / sr)
    frequencies[0] = frequencies[1]
    return fft_signal, frequencies


def xifft(fft_signal, f, multiplier=1j * 2 * np.pi):
    return np.fft.ifft(fft_signal * f * multiplier)


def dBdt_fft(B, sr=1, multiplier=1j * 2 * np.pi, sqrt=False):
    b_fft, f = xfft(B, sr)
    if sqrt:
        m = np.sqrt(multiplier * f.astype(dtype=complex))
        dBdt = np.fft.ifft(b_fft * m)
    else:
        dBdt = np.fft.ifft(b_fft * f * multiplier)
    return np.real(dBdt)


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
