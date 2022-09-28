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
import math
import shutil

from loguru import logger
import jsonschema
from jsonschema import validate
import traceback
from types import SimpleNamespace


class RecursiveNamespace(SimpleNamespace):
    """ """

    @staticmethod
    def map_entry(entry):
        if isinstance(entry, dict):
            return RecursiveNamespace(**entry)
        return entry

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        for key, val in kwargs.items():
            if type(val) == dict:
                setattr(self, key, RecursiveNamespace(**val))
            elif type(val) == list:
                setattr(self, key, list(map(self.map_entry, val)))


def frexp10(x):
    """
    Convert to mantesa exponent form
    """
    exp = int(math.log10(x))
    if exp <= 0:
        exp -= 1
    return (x / 10**exp, exp)


def frexp102str(x):
    """
    Convert to mantesa exponent form in txt
    """
    m, exp = frexp10(x)
    txt = r"%.2f$\times 10^{%d}$" % (m, exp)
    return txt


def fft(X, dT, remove_zero_frequency=True):
    """
    This function is responsible for FFT using
    numpy package of a real signal (X).
    """
    n = len(X)
    Y = 2.0 / n * np.fft.rfft(X)
    f = np.fft.rfftfreq(len(X)) / dT
    if remove_zero_frequency:
        f[0] = f[1]
    return (Y, f)


def ifft(Y):
    """
    This function is responsible for IFFT using
    numpy package of a complex FFT signal (Y).
    """
    n = len(Y)
    X = np.fft.irfft(Y) * n
    return X


def load_params(args, cfg_fname):
    """
    This method runs at initialization reads config file
    loads to argument
    """
    with open(cfg_fname, "r") as f:
        o = json.load(f)
    for k in o.keys():
        if (not hasattr(args, k)) or (vars(args)[k] is None):
            setattr(args, k, o[k])
    args.sim_id = "%03d" % args.sim_id
    args.out_dirs["base"] = args.out_dirs["base"].format(sim_id=args.sim_id)
    return args


def create_structures(args, param_file):
    """
    This method is resposible for creating
    folder structres
    """
    os.makedirs(args.out_dir, exist_ok=True)
    shutil.copy2(param_file, args.out_dir)
    return


def set_dict(x, o):
    """
    Checks if o is dictionary or not and set the associaed values
    """
    for k in o.keys():
        if isinstance(o[k], dict):
            setattr(x, k, types.SimpleNamespace())
            setattr(x, k, set_dict(getattr(x, k), o[k]))
        else:
            setattr(x, k, o[k])
    return x


def print_rec(x, spc="\t {key}->{val}"):
    """
    Print recursive attributes
    """
    for k in vars(x).keys():
        if isinstance(getattr(x, k), types.SimpleNamespace):
            print(spc.format(key=k, val=""))
            print_rec(getattr(x, k), spc="\t" + spc)
        elif isinstance(getattr(x, k), list) and isinstance(
            getattr(x, k)[0], types.SimpleNamespace
        ):
            print(spc.format(key=k, val=""))
            for ox in getattr(x, k):
                print_rec(ox, spc="\t" + spc)
                print("")
        else:
            print(spc.format(key=k, val=vars(x)[k]))
    return


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
        return "+" if sign == 0 else "-"

    files = glob.glob(base + "/*Bin*.csv")
    files.sort()
    header = (
        "* Lines starting with * are just comments.\n"
        + "* Text after the numbers is ignored \n"
        + "* BM ocean conductivity model\n"
        + "*/ %s, INF,/ ! layer thicknesses in km\n"
        + "*/ %s ,/ !Resistivities in Ohm-m\n"
        + "%d                             Number of layers from Ocean floor/Earth surface\n"
    )
    eachline = (
        "\n%.7f                      Conductivity in S/m (layer %d)"
        + "\n%.3fe%s%02d                      Layer thickness in m (layer %d)\n"
    )
    lastline = "\n1.1220100                      Semi-infinite earth conductivity"
    ocean_layer = []
    for f in files:
        o = pd.read_csv(f)
        bname = f.split("_")[-1].replace(".csv", "")
        ocean_layer.append(
            {"bin": bname, "depth": o["Thk(km)"][0] * 1e3, "rho": o["Rho(ohm-m)"][0]}
        )
        thks, rhos = ", ".join([str(x) for x in o["Thk(km)"][1:]]), ", ".join(
            [str(x) for x in o["Rho(ohm-m)"][1:]]
        )
        rhos += ", %.3f" % (1.0 / 1.1220100)
        body = header % (thks, rhos, (len(o) - 1))
        for i, row in o.iterrows():
            if i > 0:
                th = row["Thk(km)"] * 1e3
                body += eachline % (
                    1 / row["Rho(ohm-m)"],
                    i,
                    fman(th),
                    sign(th),
                    fexp(th),
                    i,
                )
        body += lastline
        with open(f.replace(".csv", ".txt"), "w") as f:
            f.writelines(body)
    ocean_layer = pd.DataFrame.from_records(ocean_layer)
    ocean_layer.to_csv(base + "/OceanLayers.csv", header=True, index=False)
    return


def validate_jsons(json_file, schema_folder, opcode):
    """
    Validate input json against the schema.
    """
    isValid, o = True, None
    with open(schema_folder + "Schema%d.json" % opcode, "r") as f:
        schema = json.load(f)
    try:
        with open(json_file, "r") as f:
            jdata = json.load(f)
        validate(instance=jdata, schema=schema)
    except Exception as e:
        traceback_str = traceback.format_exc()
        isValid = False
    if isValid:
        with open(json_file, "r") as f:
            o = json.loads(
                "\n".join(f.readlines()), object_hook=lambda d: SimpleNamespace(**d)
            )
            o.sid = "%03d" % o.sid
            o.out_dir = o.out_dir.format(sid=o.sid)
            logger.info(f"Given JSON data is Valid")
    else:
        logger.error(f"Given JSON data is InValid: {traceback_str}")
    return isValid, o

def component_mappings(field="B2E", comp="X"):
    """
    This method holds components mapping from (i) B2E
    """
    _map_ = {
        "B2E": {
            "X": "Y",
            "Y": "X"
        }
    }
    return _map_[field][comp]
    