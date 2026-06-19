import numpy as np
from loguru import logger


def sqrtf_weighted_signal(b_theta: np.ndarray, dt: float) -> np.ndarray:
    """
    Apply sqrt(f) spectral weighting to a B_theta(t) time series.

    Implements the conductivity-agnostic proxy (i2πf)^0.5 · B(f) derived in
    Section 5 (Eq. 43-44, Table 1, A=0.5 row) WITHOUT any Earth conductivity
    term — this is the proxy requested by Reviewer 1, Comment 1.1(c).

    Parameters
    ----------
    b_theta : np.ndarray
        Directional projection Bx*cos(θ) + By*sin(θ) [nT].
    dt : float
        Sample interval [s].

    Returns
    -------
    proxy : np.ndarray
        sqrt(f)-weighted time series, same length as b_theta.
    """
    N = len(b_theta)
    B_f = np.fft.rfft(b_theta)
    freqs = np.fft.rfftfreq(N, d=dt)
    weight = np.sqrt(freqs)
    weight[0] = 0.0          # zero DC component
    return np.fft.irfft(B_f * weight, n=N)


class Pipeline:
    """
    A class to represent a generic pipeline structure.
    """

    def __init__(
        self,
        angle: float = 90.0,
        ro: float = 0.381,
        ri: float = 0.3654,
        rho: float = 0.18e-6,
    ):
        """
        Initialize the pipeline with a sequence of steps.

        Parameters:
        angle (float, optional): Angle in degrees representing the orientation of the pipeline from the North angle. Default is 90 (East)
        ro (float, optional): Outer radius of the pipeline in kilometers. Default is 0.381e-3 m.
        ri (float, optional): Inner radius of the pipeline in kilometers. Default is 0.3654e-3 m.
        rho (float, optional): Resistivity of the pipeline in ohm-km. Default is 0.18e-6 ohm-m.
        """
        self.ro = ro  # Outer radius in km
        self.ri = ri  # Inner radius in km
        self.rho = rho  # Resistivity ohm-km
        self.angle = (
            angle  # Angle in degrees (oriatation of the pipeline from North angle)
        )
        self.Z = self.compute_Z(rho, ro, ri)  # Compute Z based on the given parameters
        logger.info(
            f"Pipeline initialized with angle: {self.angle} degrees, "
            f"outer radius: {self.ro} m, inner radius: {self.ri} m, "
            f"resistivity: {self.rho} ohm-m, Z: {self.Z} ohm/km"
        )
        return

    def compute_Z(self, rho: float, ro: float, ri: float) -> float:
        """
        Compute the Z value based on the resistivity and radii.

        Parameters:
        rho (float): Resistivity in ohm-km.
        ro (float): Outer radius in km.
        ri (float): Inner radius in km.

        Returns:
        float: Computed Z value in ohm/km.
        """
        Z = rho / (np.pi * (ro**2 - ri**2))  # Z in ohm/m
        Z *= 1e3  # Convert to ohm/km
        logger.debug(f"Computed Z: {Z} ohm/km using rho: {rho}, ro: {ro}, ri: {ri}")
        return Z

    def compute_E(self, ex: np.array, ey: np.array):
        """
        Compute the electric field (E-field) along the pipeline based on the given horizontal electric field components.

        Parameters:
        ex (np.array): Array of horizontal electric field values in the X-direction (North-South) in mV/km or V/km.
        ey (np.array): Array of horizontal electric field values in the Y-direction (East-West) in mV/km or V/km.

        Returns:
        np.array: Computed E-field values along the pipeline in V/km or mV/km.
        """
        logger.debug("Calculate E-field along the pipeline")
        E_pipe = (ex * np.cos(np.deg2rad(self.angle))) + (
            ey * np.sin(np.deg2rad(self.angle))
        )  # E-field in V/km or mV/km
        return E_pipe

    def compute_J(self, E: np.array = None):
        """
        Compute the current density (J/GIC) based on the electric field (E-field) values.

        Parameters:
        E (np.array): Array of electric field values in V/km or mV/km.

        Returns:
        np.array: Computed current density (J/GIC) values in A or mA.
        """
        logger.debug("Calculate GIC along the pipeline")
        gic_pipe = E if E is not None else self.E_pipe / self.Z
        return gic_pipe

    def compute_segmented_correlation(
        self,
        bx: np.array,
        by: np.array,
        ex: np.array,
        ey: np.array,
        dtheta: float = 1.0,
        trim_fraction: float = 0.3,
    ):
        """
        Compute the correlation between pipeline GIC and the magnetic field
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
        gic = self.compute_J(self.compute_E(ex, ey))
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

    def compute_proxy_robustness(
        self,
        bx: np.ndarray,
        by: np.ndarray,
        gic_A: np.ndarray,
        gic_B: np.ndarray,
        dt: float,
        dtheta: float = 1.0,
        trim_fraction: float = 0.3,
    ):
        """
        Reviewer 1, Comment 1.1(c): test whether a single sqrt(f)-weighted
        B_theta proxy — built with NO knowledge of Earth conductivity — tracks
        GIC computed under BOTH Case A and Case B across all azimuths.

        Mirrors compute_segmented_correlation() exactly (same trimming, same
        dtheta sweep).  The FFT is applied to the full (untrimmed) b_theta for
        spectral resolution, and the resulting proxy is trimmed to [i0:i1]
        before correlation, matching the trimming applied to gic_A/gic_B.

        Parameters
        ----------
        bx, by : np.ndarray
            Tapered, mean-removed Bx, By [nT] — same arrays passed to
            compute_segmented_correlation() for Figs 5/6.
        gic_A, gic_B : np.ndarray
            Full (untrimmed) GIC time series under Case A and Case B
            respectively, computed from the same bx, by input.
        dt : float
            Sample interval [s], passed to sqrtf_weighted_signal.
        dtheta, trim_fraction : float
            Match values used in compute_segmented_correlation().

        Returns
        -------
        theta_deg : np.ndarray
        r_A : np.ndarray
            |r(proxy, GIC_A)| at each azimuth.
        r_B : np.ndarray
            |r(proxy, GIC_B)| at each azimuth.
        """
        logger.debug("Computing sqrt(f) proxy robustness correlations")
        N = len(bx)
        i0 = int(N * trim_fraction)
        i1 = N - int(N * trim_fraction)
        gic_A_t = gic_A[i0:i1]
        gic_B_t = gic_B[i0:i1]

        theta_deg = np.arange(0, 360, dtheta)
        r_A = np.zeros(len(theta_deg), dtype=float)
        r_B = np.zeros(len(theta_deg), dtype=float)

        std_A = np.std(gic_A_t)
        std_B = np.std(gic_B_t)

        for i, th in enumerate(np.deg2rad(theta_deg)):
            b_theta = bx * np.cos(th) + by * np.sin(th)
            # FFT on full signal for spectral resolution, then trim
            proxy = sqrtf_weighted_signal(b_theta, dt)[i0:i1]
            if np.std(proxy) < 1e-10:
                continue
            if std_A > 1e-10:
                r_A[i] = np.abs(np.corrcoef(proxy, gic_A_t)[0, 1])
            if std_B > 1e-10:
                r_B[i] = np.abs(np.corrcoef(proxy, gic_B_t)[0, 1])

        logger.debug(
            f"Proxy robustness: peak r_A={np.max(r_A):.4f} at "
            f"{theta_deg[np.argmax(r_A)]:.1f}°, peak r_B={np.max(r_B):.4f} at "
            f"{theta_deg[np.argmax(r_B)]:.1f}°"
        )
        return theta_deg, r_A, r_B
