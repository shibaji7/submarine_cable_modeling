import numpy as np
import pandas as pd
from scipy import constants as C

from scubas.datasets import Site, PROFILES
import glob

def get_benchmark_datasets(a_scale=None):
    a_55, a_60 = 0.001*np.exp(0.115*55), 0.001*np.exp(0.115*60)
    a_scale = a_scale if a_scale else a_60/a_55
    files = glob.glob("datasets/*.csv")
    files.sort()
    print("a_scale>>>", a_scale)
    datasets = pd.concat([pd.read_csv(f, parse_dates=["datetime"]) for f in files])
    datasets.x = (datasets.x-np.mean(datasets.x.iloc[:60]))*a_scale
    datasets.y = (datasets.y-np.mean(datasets.y.iloc[:60]))*a_scale
    datasets.z = (datasets.z-np.mean(datasets.z.iloc[:60]))*a_scale
    datasets["dx"] = np.diff(datasets.x, prepend=datasets.x.iloc[0])
    datasets["dy"] = np.diff(datasets.y, prepend=datasets.y.iloc[0])
    datasets["dz"] = np.diff(datasets.z, prepend=datasets.z.iloc[0])
    return datasets

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
    gca.plot(gca.get_xlim(),gca.get_ylim(), ls="-.", lw=0.4, color="k")
    return

def abxline():
    gca = plt.gca()
    gca.set_autoscale_on(False)
    gca.plot(gca.get_xlim(),[4e-1,4e-1], ls="-.", lw=0.4, color="k")
    return


sites = [
    Site.init(
        conductivities=[1/3, 1/3000],
        thicknesses=[10000, np.inf],
        names=[
            "Sediments",
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper",
        site_name="Two layer Earth model",
    ),
    Site.init(
        conductivities=[1/3000, 1/3],
        thicknesses=[100000, np.inf],
        names=[
            "Sediments",
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper",
        site_name="Two layer Earth model",
    ),
    Site.init(
        conductivities=[1/3],
        thicknesses=[np.inf],
        names=[
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper",
        site_name="Unform Earth",
    ),
    Site.init(
        conductivities=[1/3000],
        thicknesses=[np.inf],
        names=[
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper",
        site_name="Unform Earth",
    ),
    PROFILES.DB,
]

def fft(signal, sr=1):
    fft_signal = np.fft.fft(signal)
    frequencies = np.fft.fftfreq(len(signal), 1/sr)
    frequencies[0] = frequencies[1]
    return fft_signal, frequencies

def ifft(fft_signal, f, multiplier=1j*2*np.pi):
    return np.fft.ifft(fft_signal*f*multiplier)

def dBdt_fft(B, sr=1, multiplier=1j*2*np.pi, sqrt=False):
    b_fft, f = fft(B, sr)
    if sqrt:
        m = np.sqrt(multiplier * f.astype(dtype=complex))
        dBdt = np.fft.ifft(b_fft*m)
    else:
        dBdt = np.fft.ifft(b_fft*f*multiplier)
    return np.real(dBdt)