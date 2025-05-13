"""
This method is aauming finite line current element
"""
import numpy as np
from plots import plot_real_imag_plots
from scipy import constants as C

dist = np.arange(-1000, 1001, 1) * 1e3
current_element_length = 200 * 1e3  # 100 km above ground
B = np.zeros_like(dist)
h = 100 * 1e3  # 100 km above ground
I = 1e6  # 1 Millions amps
theta_1 = 90 - (np.arctan2(dist, (current_element_length / 2)) * 180 / np.pi)
theta_2 = 90 - (np.arctan2(-dist, (current_element_length / 2)) * 180 / np.pi)
B = (
    (C.mu_0 * I)
    * (np.sin(np.deg2rad(theta_1)) + np.sin(np.deg2rad(theta_2)))
    / (4 * np.pi * h)
)
plot_real_imag_plots(dist, B, "poc_finite_current_element.png")
