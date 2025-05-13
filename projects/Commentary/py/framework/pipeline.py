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
        ex (np.array): Array of horizontal electric field values in the X-direction (North-South) in mV/km or V/km.
        ey (np.array): Array of horizontal electric field values in the Y-direction (East-West) in mV/km or V/km.
        normalize (bool, optional): Whether to normalize the correlation values. Default is False.
        rotate (float, optional): Angle in degrees to rotate the theta segments. Default is 45.0.

        Returns:
        np.array: Theta segments in degrees.
        np.array: Computed segmented correlation values.
        """
        logger.debug("Calculate Correlation Coefficients")
        theta = np.linspace(
            0, 2 * np.pi, 1001
        )  # These all are E, B, or dB/dt direction
        gic = self.compute_J(self.compute_E(ex, ey))
        cor = np.abs(
            (
                np.cos(theta - np.deg2rad(self.angle))
                + np.sin(theta - np.deg2rad(self.angle))
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
