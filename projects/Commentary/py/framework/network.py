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
        Compute the magnitude of the electric field.

        Parameters:
        ex (np.array): Array of electric field values in the X-direction (North-South).
        ey (np.array): Array of electric field values in the Y-direction (East-West).

        Returns:
        np.array: Magnitude of the electric field.
        """
        eh = np.sqrt(ex**2 + ey**2)
        return eh

    def compute_J(
        self,
        ex: Sequence[float],
        ey: Sequence[float],
    ):
        """
        Calculate the geomagnetically induced current (GIC) at this substation
        based on the electric field components.

        GIC = a * Ex + b * Ey

        where a, b are the substation coupling coefficients from the network
        admittance model.

        Parameters:
        ex (Sequence[float]): Electric field X-component (North-South).
        ey (Sequence[float]): Electric field Y-component (East-West).

        Returns:
        np.array: GIC values in amperes.
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
        bx: np.array,
        by: np.array,
        ex: np.array,
        ey: np.array,
        dtheta: float = 1.0,
        trim_fraction: float = 0.1,
    ):
        """
        Compute the correlation between substation GIC and the magnetic field
        projection B_theta(t) = Bx(t)*cos(theta) + By(t)*sin(theta) as a
        function of azimuthal direction theta.

        This follows the approach of Trichtchenko and Boteler (2004), where
        the directional dependence of correlation is determined by projecting
        the measured Bx and By components onto each azimuthal direction and
        computing the Pearson correlation with GIC at each angle.

        Parameters:
        bx (np.array): Array of north (X) magnetic field (or dB/dt) variations.
        by (np.array): Array of east (Y) magnetic field (or dB/dt) variations.
        ex (np.array): Array of north (X) electric field values in mV/km.
        ey (np.array): Array of east (Y) electric field values in mV/km.
        dtheta (float): Angular resolution in degrees. Default is 1.0.
        trim_fraction (float): Fraction of time series to trim from each end to
            exclude tapered/quiet boundary regions. Default is 0.1 (10% each end).

        Returns:
        theta_deg (np.array): Azimuthal angles in degrees [0, 360).
        cor (np.array): Correlation coefficient r(GIC, B_theta) at each angle.
        """
        logger.debug("Calculate directional correlation coefficients")

        # Trim tapered boundary regions before computing correlations.
        # The benchmark B-field has a 5% cosine taper on each end; the gradient
        # (dB/dt) has large artificial transients there that distort the result.
        N = len(bx)
        i_start = int(N * trim_fraction)
        i_end = N - int(N * trim_fraction)
        bx_t = bx[i_start:i_end]
        by_t = by[i_start:i_end]

        # Compute GIC over full arrays, then trim to the same active window
        gic = self.compute_J(ex, ey)
        gic_t = gic[i_start:i_end]

        theta_deg = np.arange(0, 360, dtheta)
        theta_rad = np.deg2rad(theta_deg)

        cor = np.zeros(len(theta_deg))
        for i, th in enumerate(theta_rad):
            # Project B-field (or dB/dt) onto direction theta
            b_theta = bx_t * np.cos(th) + by_t * np.sin(th)
            if np.std(b_theta) < 1e-10 or np.std(gic_t) < 1e-10:
                cor[i] = 0.0
            else:
                cor[i] = np.corrcoef(b_theta, gic_t)[0, 1]

        logger.debug(
            f"Peak |r| = {np.max(np.abs(cor)):.4f} at "
            f"theta = {theta_deg[np.argmax(np.abs(cor))]:.1f} deg"
        )
        return theta_deg, cor
