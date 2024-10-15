"""
    info_theory.py: 
        Module has multiple functions to compute information and mutual
        information between two datasets (discreat and analog). Most of the
        parameters are taken from the existing studies. Following refrences 
        provide a detailed description of the parameters and give a brief 
        description of how to numerically compute them.

    * Wing, S., Johnson, J. R., Turner, D. L., Ukhorskiy, A. Y., & Boyd, A. J. (2022). 
        Untangling the solar wind and magnetospheric drivers of the radiation belt electrons.
        Journal of Geophysical Research: Space Physics, 127, e2021JA030246. 
        https://doi.org/10.1029/2021JA030246
    * Wing, S., J. R. Johnson, E. Camporeale, and G. D. Reeves (2016), Information theoretical
        approach to discovering solar wind drivers of the outer radiation belt, 
        J. Geophys. Res. Space Physics, 121, 9378–9399, doi:10.1002/2016JA022711.
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S.; Coyle, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde


class analog(object):
    """
    This class holds static methods to compute information theory based parameters
    for analog signals.
    """

    @staticmethod
    def differential_entropy(x: np.array, N: int = 100) -> float:
        """
        Calculate the differential entropy of a continuous random variable.
            H(x), equation (2) in citation.

        Attributes:
        -----------
            x (array-like): Input data sample
            N (int): Number of interpolated values

        Returns:
        --------
            S (float): Differential entropy
        """
        # Use Gaussian kernel density estimation to estimate the probability
        # density function
        kde = gaussian_kde(x)
        # Evaluate the PDF on a range of values
        values = np.linspace(min(x), max(x), N)
        pdf = kde(values)
        # Calculate the differential entropy
        S = -np.sum(pdf * np.log(pdf) * (values[1] - values[0]))
        return S

    @staticmethod
    def joint_differential_entropy(x: np.array, y: np.array, N: int = 100):
        """
        Calculate the joint differential entropy of two continuous random variables.
            H(x,y), equation (3) in citation.

        Attributes:
        -----------
            x (array-like): Input data sample
            y (array-like): Input data sample
            N (int): Number of interpolated values

        Returns:
        --------
            S (float): Joint differential entropy
        """
        # Stack the data for joint KDE
        data = np.vstack([x, y])
        # Use Gaussian kernel density estimation to estimate the joint PDF
        kde = gaussian_kde(data)
        # Create a grid of points to evaluate the joint PDF
        xx = np.linspace(min(x), max(x), N)
        yy = np.linspace(min(y), max(y), N)
        X, Y = np.meshgrid(xx, yy)
        positions = np.vstack([X.ravel(), Y.ravel()])
        # Evaluate the joint PDF on the grid
        joint_pdf = kde(positions).reshape(X.shape)
        # Calculate the joint differential entropy
        dx = xx[1] - xx[0]
        dy = yy[1] - yy[0]
        S = -np.sum(joint_pdf * np.log(joint_pdf) * dx * dy)
        return S

    @staticmethod
    def conditional_differential_entropy(x: np.array, y: np.array, N: int = 100):
        """
        Calculate the conditional differential entropy H(x|y).
            H(x|y), equations in citation.

        Attributes:
        -----------
            x (array-like): Input data sample
            y (array-like): Input data sample
            N (int): Number of interpolated values

        Returns:
        --------
            H_x_y (float): Conditional differential entropy
        """
        # Joint differential entropy
        H_xy = analog.joint_differential_entropy(x, y, N=N)
        # Differential entropy of y
        H_y = analog.differential_entropy(y)
        # Conditional differential entropy
        H_x_y = H_xy - H_y
        return H_x_y

    @staticmethod
    def mutual_information(x: np.array, y: np.array, N: int = 100):
        """
        Calculate the mutual information between two continuous random variables.
            MI(x,y) = H(x) + H(y) - H(x,y)

        Attributes:
        -----------
            x (array-like): Input data sample
            y (array-like): Input data sample
            N (int): Number of interpolated values

        Returns:
        --------
            MI (float): Mutual information between x, y.
        """
        # Calculate entropy of x, H(x)
        H_x = analog.differential_entropy(x, N=N)
        # Calculate entropy of y, H(y)
        H_y = analog.differential_entropy(y, N=N)
        # Calculate joint entropy of x and y, H(x, y)
        H_xy = analog.joint_differential_entropy(x, y, N=N)
        # Calculate mutual information
        mi = H_x + H_y - H_xy
        return mi

    @staticmethod
    def conditional_mutual_information(
        x: np.array, y: np.array, Z: np.array, N: int = 100
    ):
        """
        Calculate the conditional mutual information CI(x,y|Z).
            CI(x,y|Z) = H(x|Z) + H(y|Z) - H(x,y|Z)

        Attributes:
        -----------
            x (array-like): Input data sample
            y (array-like): Input data sample
            Z (array-like): Input data sample (condition upon)
            N (int): Number of interpolated values

        Returns:
        --------
            CI (float): Conditional mutual information between x, y given Z.
        """
        # Calculate conditional entropy H(x|Z)
        H_x_Z = analog.conditional_differential_entropy(x, Z, N=N)
        # Calculate conditional entropy H(y|Z)
        H_y_Z = analog.conditional_differential_entropy(y, Z, N=N)
        # Calculate conditional entropy H(x,y|Z)
        H_xy_Z = analog.conditional_differential_entropy(np.vstack([x, y]), Z, N=N)
        # Calculate conditional information H(x|Z) + H(y|Z) - H(x,y|Z)
        CI = H_x_Z + H_y_Z - H_xy_Z
        return CI

    @staticmethod
    def transfer_entropy(x: np.array, y: np.array, N: int = 100, delay: int = 1):
        """
        Calculate the transfer entropy from data_source (x) to data_target (y) TE(x->y).
            TE(x->y) = CMI(y[t+d],x[t]|y_p[t])

        Attributes:
        -----------
            x (array-like): Input data sample
            y (array-like): Input data sample
            N (int): Number of interpolated values
            delay (int): delay to check the timeseries

        Returns:
        --------
            TE (float): A special case conditional mutual information between x and y.
        """
        x_past, y_past = (x[:-delay], y[:-delay])
        y_future = y[delay:]
        # Calculate H(y[t+d]|y_p[t])
        H_yfuture_ypast = analog.conditional_differential_entropy(y_future, y_past)
        # Calculate H(y[t+d]|y_p[t],x[t])
        H_yfuture_ypast_xpast = analog.conditional_differential_entropy(
            y_future, np.vstack([y_past, x_past])
        )
        # Calculate TE(x->y) = CMI(y[t+d],x[t]|y_p[t]) = H(y[t+d]|y_p[t]) - H(y[t+d]|y_p[t],x[t])
        TE = H_yfuture_ypast - H_yfuture_ypast_xpast
        return TE


class discrete(object):
    """
    This class holds static methods to compute information theory based parameters
    for discrete signals.
    """

    @staticmethod
    def entropy(
        pk: np.array,
        qk: np.array = None,
        base: int = 2,
        axis: int = 0,
        nan_policy: str = "propagate",
        keepdims: bool = False,
    ):
        """
        Calculate the Shannon entropy/relative entropy of given distribution(s).
            * If only probabilities pk are given, the Shannon entropy is calculated as
                H = -sum(pk * log(pk)).
            * If qk is not None, then compute the relative entropy D = sum(pk * log(pk / qk)).
                This quantity is also known as the Kullback-Leibler divergence.
            * This routine will normalize pk and qk if they don’t sum to 1.

        Attributes:
        -----------
            pk (array_like): Defines the (discrete) distribution. Along each axis-slice of pk,
                element i is the (possibly unnormalized) probability of event i.
            qk (array_like, opt): Sequence against which the relative entropy is computed.
                Should be in the same format as pk.
            base (float, opt): The logarithmic base to use, defaults to e (natural logarithm).
            axis (int, None): The axis of the input along which to compute the statistic.
            nan_policy: Check https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.htm.
            keepdims: Check https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.htm.

        Returns:
        --------
            S (float, array_like): The calculated entropy.
        """
        return stats.entropy(
            pk, qk=qk, base=base, axis=axis, nan_policy=nan_policy, keepdims=keepdims
        )

    @staticmethod
    def joint_entropy(
        p_xy: np.array,
        base: int = 2,
    ):
        """
        Calculate the joint entropy H(X, Y) given the joint probability distribution p(x, y).

        Attibutes:
        ----------
            p_xy (2D array-like): Joint probability distribution p(x, y).
            base (float): Log base to calculate

        Returns:
            S (float): Joint entropy H(X, Y).
        """
        p_xy = np.array(p_xy)
        # Remove zero probabilities to avoid log(0)
        p_xy = p_xy[p_xy >= 0]
        S = -np.nansum(p_xy * np.emath.logn(base, p_xy))
        return S

    @staticmethod
    def conditional_entropy(
        p_xy: np.array,
        p_x: np.array = None,
        base: int = 2,
    ):
        """
        Calculate the conditional entropy H(Y|X).

        Attibutes:
        ----------
            p_xy (2D array-like): Joint probability distribution p(x, y).
            p_x (array-like): 1D array of marginal probabilities p(x)
            base (float): Log base to calculate

        Returns:
            S (float): Conditional entropy H(Y|X).
        """
        # If p_x is not given estimate from p(x,y)
        p_x = np.sum(p_xy, axis=1) if p_x is None else p_x
        H_y_x = 0.0
        for i in range(p_xy.shape[0]):  # Loop through p(x)
            for j in range(p_xy.shape[1]):  # Loop through p(y)
                if p_xy[i, j] > 0 and p_x[i] > 0:  # Integrate through p(x), p(x,y) > 0
                    H_y_x -= p_xy[i, j] * np.emath.logn(base, p_xy[i, j] / p_x[i])
        return H_y_x

    @staticmethod
    def mutual_information(
        p_xy: np.array,
        p_x: np.array = None,
        p_y: np.array = None,
        base: int = 2,
    ):
        """
        Calculate the mutual information MI(x,y)

        Attibutes:
        ----------
            p_xy (2D array-like): Joint probability distribution p(x, y).
            p_x (array-like): 1D array of marginal probabilities p(x)
            p_y (array-like): 1D array of marginal probabilities p(y)
            base (float): Log base to calculate

        Returns:
            mi (float): Mutual Information.
        """
        # If p_x is not given estimate from p(x,y)
        p_x = np.sum(p_xy, axis=1) if p_x is None else p_x
        # If p_y is not given estimate from p(x,y)
        p_y = np.sum(p_xy, axis=0) if p_y is None else p_y
        mi = 0.0
        for i in range(p_xy.shape[0]):  # Loop through p(x)
            for j in range(p_xy.shape[1]):  # Loop through p(y)
                if p_xy[i, j] > 0:  # Integrate through p(x,y) > 0
                    mi += p_xy[i, j] * np.emath.logn(
                        base, p_xy[i, j] / (p_x[i] * p_y[j])
                    )
        return mi

    @staticmethod
    def conditional_mutual_information(
        p_xyz: np.array,
        p_xz: np.array = None,
        p_yz: np.array = None,
        p_z: np.array = None,
        base: int = 2,
    ):
        """
        Calculate the conditional mutual information CI(x,y|Z).
            CI(x,y|Z) = H(x|Z) + H(y|Z) - H(x,y|Z)

        Attibutes:
        ----------
            p_xyz (3D array-like): Joint probability distribution p(x, y, z).
            p_xz (2D array-like): 2D array of marginal probabilities p(x, z)
            p_yz (2D array-like): 2D array of marginal probabilities p(y, z)
            p_z (array-like): 1D array of marginal probabilities p(z)
            base (float): Log base to calculate

        Returns:
            cmi (float): Conditional Mutual Information.
        """
        # If p_xz is not given estimate from p(x,y,z)
        p_xz = np.sum(p_xyz, axis=1) if p_xz is None else p_xz
        # If p_yz is not given estimate from p(x,y,z)
        p_yz = np.sum(p_xyz, axis=0) if p_yz is None else p_yz
        # If p_z is not given estimate from p(x,y,z)
        p_z = np.sum(p_xyz, axis=(0, 1)) if p_z is None else p_z
        cmi = 0.0
        for i in range(p_xyz.shape[0]):  # Loop through p(x)
            for j in range(p_xyz.shape[1]):  # Loop through p(y)
                for k in range(p_xyz.shape[2]):  # Loop through p(z)
                    if (
                        p_xyz[i, j, k] > 0 and p_z[k] > 0
                    ):  # Integrate through p(x,y,z), p(z) > 0
                        pr_xyz = p_xyz[i, j, k]
                        pr_xz = p_xz[i, k]
                        pr_yz = p_yz[j, k]
                        pr_z = p_z[k]
                        pr_x_z = pr_xz / pr_z
                        pr_y_z = pr_yz / pr_z
                        pr_xy_z = pr_xyz / pr_z
                        cmi += p_xyz * np.emath.logn(base, pr_xy_z / (pr_x_z * pr_y_z))
        return cmi

    @staticmethod
    def compute_probabilities(x, y):
        return

    @staticmethod
    def quantize_continuous_signal(x):
        return
