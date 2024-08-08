import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import numpy as np
import pandas as pd
from scipy import constants as C

f = 1e-3
rho1, rho2 = 3, 3000
sig1, sig2 = 1/rho1, 1/rho2
l1 = 10000
k1 = np.sqrt(1j*2*np.pi*f*C.mu_0*sig1)
eta1 = 1j*2*np.pi*f/k1
K2 = 1j*2*np.pi*f/np.sqrt(1j*2*np.pi*f*C.mu_0*sig2)
r2 = (eta1 - K2) / (eta1 + K2)

print(r2)
