import numpy as np
import pandas as pd
from intermagnet import Efield
from loguru import logger
from scipy import constants as C
from scipy.signal import get_window
from scubas.datasets import PROFILES, Site
from scubas.utils import fft, ifft

setattr(
    PROFILES,
    "CaseA",
    Site.init(
        conductivities=[1 / 3, 1 / 3000],
        thicknesses=[10000, np.inf],
        names=[
            "Sediments",
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper",
        site_name="Two layer Earth model",
    ),
)
setattr(
    PROFILES,
    "CaseB",
    Site.init(
        conductivities=[1 / 3000, 1 / 3],
        thicknesses=[100000, np.inf],
        names=[
            "Sediments",
            "Lower Mantle",
        ],
        desciption="This model is modified for David's cometary paper",
        site_name="Two layer Earth model",
    ),
)
setattr(
    PROFILES,
    "Uniform",
    Site.init(
        conductivities=[1 / 3],
        thicknesses=[np.inf],
        names=[
            "Sediments",
        ],
        desciption="This model is modified for David's cometary paper",
        site_name="1 layer Uniform Earth model",
    ),
)


class TransferFunction:
    """
    A class to represent a transfer function.
    """

    def __init__(
        self,
        site: Site,
    ):
        """
        Initialize the transfer function with frequency, site, and data.
        :param site: Site object containing the model
        """
        self.site = site
        return

    def get_tapering_function(self, t, p):
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

    def calcZ(self):
        omega = 2 * np.pi * self.freqs
        sigma_s = 1 / self.site.layers[0].resistivity
        k2 = 1.0j * omega * C.mu_0 * sigma_s
        k = np.sqrt(k2)
        Zo = (1.0j * omega * C.mu_0 / k) / (C.mu_0 / 1.0e-3)
        return Zo

    def calcTF(self, frequencies: np.array = None):
        """
        Calculate the transfer functions.
        """
        Zo = self.calcZ()
        Zd = self.site.calcZ(self.freqs)[1]
        omega = 2 * C.pi * self.freqs
        sigma_s = 1 / self.site.layers[0].resistivity
        k2 = 1.0j * omega * C.mu_0 * sigma_s
        k = np.sqrt(k2)
        kd = k * 0

        func = Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Zo))
        tf = pd.DataFrame()
        tf["freq"], tf["E2B"] = self.freqs, func
        if frequencies is not None:
            for f in frequencies:
                findex = np.abs(self.freqs - f).argmin()
                logger.debug(
                    f"Frequency: {self.freqs[findex]} -> {np.angle(tf.E2B[findex], deg=True).__format__('.2f')} / {np.abs(tf.E2B[findex]).__format__('.4f')}"
                )
        return tf

    def compute_Efield(self, bx: np.array, by: np.array, del_t: float, p: float = 0.05):
        """
        Compute the electric field.
        :param data: Data object containing the electric field data
        """
        window = self.get_tapering_function(np.arange(len(bx)), p=p)
        bxf, self.freqs = fft(bx * window, del_t)
        byf, _ = fft(by * window, del_t)
        t = self.calcTF()
        ey, ex = -ifft(np.array(t.E2B) * bxf), ifft(np.array(t.E2B) * byf)
        ef = Efield(ex=ex, ey=ey)
        return ef

    def find_efrq(
        self, bx: np.array, by: np.array, del_t: float, frequencies: np.array
    ):
        """
        Find the electric field in freq domain.
        :param data: Data object containing the electric field data
        :return: Electric field in frequency domain
        """
        window = get_window("hann", len(bx))
        bxf, self.freqs = fft(bx * window, del_t)
        byf, _ = fft(by * window, del_t)
        t = self.calcTF(frequencies)
        eyf, exf = np.array(t.E2B) * bxf, np.array(t.E2B) * byf
        for f in frequencies:
            findex = np.abs(self.freqs - f).argmin()
            logger.debug(
                f"Frequency: {self.freqs[findex]} -> eyf[Index]: {np.angle(eyf[findex], deg=True).__format__('.2f')} / {np.abs(eyf[findex]).__format__('.2f')}"
            )
        for f in frequencies:
            findex = np.abs(self.freqs - f).argmin()
            logger.debug(
                f"Frequency: {self.freqs[findex]} -> bxf[Index]: {np.angle(bxf[findex], deg=True).__format__('.2f')} / {np.abs(bxf[findex]).__format__('.2f')}"
            )
        return
