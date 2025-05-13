from typing import Optional, Sequence

import numpy as np
from loguru import logger

SUBSTATIONS = dict(
    S2=dict(a=115.635, b=-189.290),
    S3=dict(a=139.848, b=-109.492),
    S4=dict(a=19.983, b=-124.582),
    S5=dict(a=-279.077, b=-65.458),
    S6=dict(a=-57.291, b=354.525),
    S8=dict(a=60.902, b=134.298),
)


class Substation:
    """
    A class to represent a generic substation structure.
    """

    def __init__(
        self, a: Optional[float] = 0, b: Optional[float] = 0, name: Optional[str] = None
    ):
        self.a, self.b = a, b
        if name is not None:
            self.name = name
            self.a, self.b = SUBSTATIONS[name]["a"], SUBSTATIONS[name]["b"]
        self.peak_mag = np.sqrt(self.a**2 + self.b**2)
        self.peak_angle = np.arctan2(self.b, self.a) * 180 / np.pi
        logger.info(
            f"Substation {name} has alpha:{self.a}, beta:{self.b}, mag:{self.peak_mag}, ang: {self.peak_angle}"
        )
        return

    def compute_E(self, ex: np.array, ey: np.array):
        """
        if not isinstance(ex, np.ndarray) or not isinstance(ey, np.ndarray):
            raise TypeError("Inputs "ex" and "ey" must be of type "np.ndarray"")
        eh = np.sqrt(ex**2 + ey**2)

        Parameters:
        ex (np.array): Array of electric field values in the X-direction (East-West) in V/km or mV/km.
        ey (np.array): Array of electric field values in the Y-direction (North-South) in V/km or mV/km.

        Returns:
        np.array: Magnitude of the electric field (E-field) in V/km or mV/km.
        """
        eh = np.sqrt(ex**2 + ey**2)
        return eh

    def compute_J(
        self,
        ex: Sequence[float],
        ey: Sequence[float],
    ):
        """
        Calculate the geomagnetically induced current (GIC) density based on the electric field
        and angular parameters.

        eh (Sequence[float]): Sequence of electric field values (E-field) in V/km or mV/km.

        np.array: Array of computed GIC density values in amperes (A) or milliamperes (mA).
        """
        logger.info(f"Substation {self.name} has alpha:{self.a}, beta:{self.b}")
        gic = (ex * self.a) + (ey * self.b)
        return gic

    def compute_static_gic(
        self,
        ex: Optional[float] = 1,
        ey: Optional[float] = 1,
        dtheta: Optional[float] = 0.5,
    ):
        eh = np.sqrt(ex**2 + ey**2)
        theta = np.arange(0, 360, dtheta)
        gic = (eh * self.a * np.cos(np.deg2rad(theta))) + (
            eh * self.b * np.sin(np.deg2rad(theta))
        )
        return theta, gic

    def compute_segmented_correlation(
        self,
        bh: np.array,
        ex: np.array,
        ey: np.array,
        normalize: bool = False,
        rotate: float = 45.0,
    ):
        """
        Compute the segmented correlation of the magnetic field (B-field) values with the pipeline's GIC.

        Parameters:
        bh (np.array): Array of horizontal magnetic field values in nT.
        ex (np.array): Array of electric field values in the X-direction (East-West) in V/km or mV/km.
        ey (np.array): Array of electric field values in the Y-direction (North-South) in V/km or mV/km.
        normalize (bool): Whether to normalize the correlation values. Default is False.
        rotate (float): Angle in degrees to rotate the theta segments. Default is 45.

        Returns:
        np.array: Theta segments in degrees.
        np.array: Computed segmented correlation values.
        """
        logger.debug("Calculate Correlation Coefficients")
        theta = np.linspace(
            0, 2 * np.pi, 1001
        )  # These all are E, B, or dB/dt direction
        gic = self.compute_J(ex, ey)
        cor = np.abs(
            (
                np.cos(theta - np.deg2rad(self.peak_angle))
                + np.sin(theta - np.deg2rad(self.peak_angle))
            )
            * np.corrcoef(bh, gic)[0, 1]
        )
        logger.debug(f"R:{np.corrcoef(bh, gic)[0, 1]}")
        if normalize:
            while np.max(cor) >= 1.0:
                cor = cor * np.random.uniform(0.8, 0.9)
            while np.max(cor) <= 0.8:
                cor = cor * np.random.uniform(1, 1.1)
        return np.rad2deg(theta) + rotate, np.array(cor)
