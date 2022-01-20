"""utility.py: Module is used to implement utility methods"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import numpy as np
import json
import os

def create_synthetic_B_field(Am, Phim, Tm, t=None):
    """
    This function is responsible for creating a 
    synthetic magnetic (B) field with following sets 
    of parameters.
    
    Parameter:
    ----------
    Am (list) - Magntitue at different freuency components (m)
    Phim (list) - Phase at different freuency components (m)
    Tm (list) - Periods of different freuency components (m)
    t (list) - Total time
    """
    t = np.linspace(0,60*60*72,60*60*72,endpoint=False) if t == None else t
    Bt = np.zeros(len(t))
    for A, Phi, T in zip(Am, Phim, Tm):
        Bt += A*np.sin(2*np.pi*t/T + np.deg2rad(Phi))
    return Bt, t

def fft(X, dT, remove_zero_frequency=True):
    """
    This function is responsible for FFT using 
    numpy package of a real signal (X).
    """
    n = len(X)
    Y = 2.0/n * np.fft.rfft(X)
    f = np.fft.rfftfreq(len(X))/dT
    if remove_zero_frequency: f[0] = f[1]
    return (Y, f)

def ifft(Y):
    """
    This function is responsible for IFFT using 
    numpy package of a complex FFT signal (Y).
    """
    n = len(Y)
    X = np.fft.irfft(Y)*n/2
    return X

def load_params(args, cfg_fname="py/config/parameters.json"):
    """
    This method runs at initialization reads config file
    loads to argument
    """
    with open(cfg_fname, "r") as f: o = json.load(f)
    for k in o.keys():
        if (not hasattr(args, k)) or (vars(args)[k] is None): setattr(args, k, o[k])
    args.sim_id = "%03d"%args.sim_id
    args.out_dirs["base"] = args.out_dirs["base"].format(sim_id=args.sim_id)
    args.out_dirs["synthetic"] = args.out_dirs["base"] + args.out_dirs["synthetic"]
    return args

def create_structures(args):
    """
    This method is resposible for creating 
    folder structres
    """
    os.makedirs(args.out_dirs["base"], exist_ok=True)
    os.makedirs(args.out_dirs["synthetic"], exist_ok=True)
    os.makedirs("/".join(args.in_file.split("/")[:-1]), exist_ok=True)
    return

def set_dict(x, o):
    """
    Checks if o is dictionary or not and set the associaed values
    """
    for k in o.keys():
        if isinstance(o[k], dict): x = set_dict(x, o[k])
        else: setattr(x, k, o[k])
    return x

def get_tapering_function(t, p=0.1):
    """
    This method is resposible for generateing 
    tapering function based on time sequence t 
    and tapering coefficient p
    """
    T = len(t)
    P, P2 = int(T*p), int(T*p/2)
    w = np.zeros_like(t)
    w[:P2] = 0.5*(1 - np.cos(2*np.pi*t[:P2]/P))
    w[P2:T-P2] = 1.
    w[T-P2:] = 0.5*(1 - np.cos(2*np.pi*(t[-1]-t[T-P2:])/P))
    return w