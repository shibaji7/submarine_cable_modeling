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
import glob
from decimal import Decimal
import pandas as pd
import types

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
        if isinstance(o[k], dict): 
            setattr(x, k, types.SimpleNamespace())
            setattr(x, k, set_dict(getattr(x, k), o[k]))
        else: setattr(x, k, o[k])
    return x

def print_rec(x, spc="\t {key}->{val}"):
    """
    Print recursive attributes
    """
    for k in vars(x).keys():
        if isinstance(getattr(x, k), types.SimpleNamespace): 
            print(spc.format(key=k, val=""))
            print_rec(getattr(x, k), spc="\t"+spc)
        else: print(spc.format(key=k, val=vars(x)[k]))
    return

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

def detrend_magnetic_field(B, t, p=0.1):
    """
    This method is resposible for detrend
    magnetic field data and taper it to reduce
    spurious frequency components.
    """
    w = get_tapering_function(t, p)
    B = B*w
    return B

def toBEZpy(base="data/OceanModels/"):
    """
    This method is dedicated to convert the csv file to 
    BEZpy readable text files. All the .csv files under this 
    base location will be converted. 
    """
    def fexp(number):
        (sign, digits, exponent) = Decimal(number).as_tuple()
        return len(digits) + exponent - 1
    
    def fman(number):
        return Decimal(number).scaleb(-fexp(number)).normalize()
    
    def sign(number):
        (sign, digits, exponent) = Decimal(number).as_tuple()
        return "+" if sign==0 else "-"
    
    files = glob.glob(base+"/*Bin*.csv")
    files.sort()
    header = "* Lines starting with * are just comments.\n"+\
                "* Text after the numbers is ignored \n"+\
                "* BM ocean conductivity model\n"+\
                "*/ %s, INF,/ ! layer thicknesses in km\n"+\
                "*/ %s ,/ !Resistivities in Ohm-m\n"+\
                "%d                             Number of layers from surface\n"
    eachline = "\n%.7f                      Conductivity in S/m (layer %d)"+\
                "\n%.3fe%s%02d                      Layer thickness in m (layer %d)\n"
    lastline = "\n1.1220100                      Semi-infinite earth conductivity"
    ocean_layer = []
    for f in files:
        o = pd.read_csv(f)
        bname = f.split("_")[-1].replace(".csv","")
        ocean_layer.append({"bin": bname, "depth": o["Thk(km)"][0]*1e3, "rho": o["Rho(ohm-m)"][0]})
        thks, rhos = ", ".join([str(x) for x in o["Thk(km)"][1:]]),\
                ", ".join([str(x) for x in o["Rho(ohm-m)"][1:]])
        rhos += (", %.3f"%(1./1.1220100))
        body = header%(thks, rhos, (len(o)-1))
        for i, row in o.iterrows():
            if i > 0:
                th = row["Thk(km)"]*1e3
                body += eachline%(1/row["Rho(ohm-m)"], i, fman(th), sign(th), fexp(th), i)
        body += lastline
        with open(f.replace(".csv", ".txt"), "w") as f: f.writelines(body)
    ocean_layer = pd.DataFrame.from_records(ocean_layer)
    ocean_layer.to_csv(base + "/OceanLayers.csv", header=True, index=False)
    return