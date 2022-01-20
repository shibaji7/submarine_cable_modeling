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

import utility
from OM import SynB

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
        synb = SynB(self.B_syn, self.tp, self.verbose, 
                    self.out_dirs["synthetic"])
        synb.run()
        return
    
    @staticmethod
    def run(args):
        """
        This is a static method that helps to run
        the simulation step-by-step manner
        """
        sim = Simulation(args)
        if sim.syn: sim.synthetic_B_field_simulation()
        return
            

# Script run can also be done via main program
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-sid", "--sim_id", default=None, help="Simulation ID, default in params.json")
    parser.add_argument("-sy", "--syn", action="store_true", help="Run synthetic B-field analysis")
    parser.add_argument("-v", "--verbose", action="store_false", help="Increase output verbosity (default True)")
    args = parser.parse_args()
    logger.info(f"Simulation run using fitacf_amgeo.__main__")
    if args.verbose:
        logger.info("Parameter list for simulation ")
        for k in vars(args).keys():
            print("     ", k, "->", vars(args)[k])
    logger.info(f"Loading parameters from cfg_prm.json")
    args = utility.load_params(args)
    logger.info(f"Fork folder structure to save data/figs.")
    utility.create_structures(args)
    Simulation.run(args)
    logger.info(f"Simulation end! Clear local files.")
    os.system("rm -rf `find -type d -name .ipynb_checkpoints`:")
    os.system("rm -rf __pycache__/")
    pass