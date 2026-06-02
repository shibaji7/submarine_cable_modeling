import numpy as np
from loguru import logger


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
