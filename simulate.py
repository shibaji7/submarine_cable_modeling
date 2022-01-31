"""simulate.py: Module is used to implement all simulations avalable under this codebase"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
import sys
sys.path.extend(["py/", "py/config/"])
import datetime as dt
import argparse
from dateutil import parser as prs
from loguru import logger
import json

import utility
from oml import SynB
from cable import EventAnalysis

class Simulation(object):
    """
    Class is dedicated for all simulation runs
    """
    
    def __init__(self, args):
        for p in vars(args).keys():
            setattr(self, p, vars(args)[p])
        return
    
    def synthetic_B_field_simulation(self):
        """
        Create and invoke synthetic B-field simulation only
        """
        synb = SynB(self.B_syn, self.tp, self.verbose, self.out_dirs["synthetic"])
        synb.run()
        return
    
    def analyze_event_data(self):
        """
        Create and invoke event based analysis.
        Stucture of the code follows:
        a. Any cable can be sub-devided multiple segment along the path
        b. Data and Ocean models will be provided in data/OceanModels/ dirs
        c. Calculate TF and E(t) for floor and surface.
        d. Calculate induced V(t) along the cable section and total voltage (sum).
        e. DSTL voltage calculation
        f. Total voltage calculation
        
        Preprocessings:
        ---------------
        1. Read mapping file
        2. To Ben's csv file to BEZpy reaable .txt
        """
        logger.info(f"Load mapping {self.in_file}")
        with open(self.in_file, "r") as f: o = json.load(f)
        logger.info(f"Convert to BEZpy")
        utility.toBEZpy("/".join(o["model_location"].split("/")[:-1]))
        logger.info(f"Start simulation....")
        e = EventAnalysis(o, self.out_dirs["base"], self.verbose, 
                          self.syne, self.E_syn, self.out_dirs["synthetic"])
        e.calclulate_total_parameters()
        logger.info(f"Simulation end!")
        return
    
    @staticmethod
    def run(args):
        """
        This is a static method that helps to run
        the simulation step-by-step manner
        """
        sim = Simulation(args)
        if sim.synb: sim.synthetic_B_field_simulation()
        sim.analyze_event_data()
        return
            

# Script run can also be done via main program
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-sid", "--sim_id", default=None, help="Simulation ID, default in params.json")
    parser.add_argument("-syb", "--synb", action="store_true", help="Run synthetic B-field analysis")
    parser.add_argument("-sye", "--syne", action="store_true", help="Run synthetic E-field analysis")
    parser.add_argument("-v", "--verbose", action="store_false", help="Increase output verbosity (default True)")
    args = parser.parse_args()
    logger.info(f"Simulation run using fitacf_amgeo.__main__")
    if args.verbose:
        logger.info("Parameter list for simulation ")
        utility.print_rec(args)
    logger.info(f"Loading parameters from cfg_prm.json")
    args = utility.load_params(args)
    logger.info(f"Fork folder structure to save data/figs.")
    utility.create_structures(args)
    Simulation.run(args)
    logger.info(f"Simulation end! Clear local files.")
    os.system("rm -rf `find -type d -name .ipynb_checkpoints`:")
    os.system("rm -rf __pycache__/")
    pass
