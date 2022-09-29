"""simulate.py: Module is used to implement all simulations avalable under this codebase"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import argparse
import os

from loguru import logger

from gic.model import utils
from gic.validation import SytheticCableAnalysis as SCA

# from synthetic import MagneticFieldAnalysis as SynB
# #from cables import Cable


class Simulation(object):
    """
    Class is dedicated for all simulation runs
    """

    def __init__(self, args):
        self.args = args
        return

    def B_field_simulation(self):
        """
        Create and invoke synthetic B-field simulation only
        """
        synb = SCA(self.args.Bfield, self.args.cable, out_dir=self.args.out_dir)
        synb.run()
        #         if "cable" in self.args.__dict__.keys():
        #             cab = Cable(self.args, None, synb.Bfield, synb.components)
        #             cab.run_nodal_analysis()
        return

    #     def E_field_simulation(self):
    #         """
    #         Create and invoke synthetic E-field simulation only
    #         """
    #         syne = SynE(self.args)
    #         syne.run()
    #         if "cable" in self.args.__dict__.keys():
    #             cab = Cable(self.args, syne.Efield, None, syne.components)
    #             cab.run_nodal_analysis()
    #         return

    #     #TODO
    #     def analyze_event_data(self):
    #         """
    #         Create and invoke event based analysis.
    #         Stucture of the code follows:
    #         a. Any cable can be sub-devided multiple segment along the path
    #         b. Data and Ocean models will be provided in data/OceanModels/ dirs
    #         c. Calculate TF and E(t) for floor and surface.
    #         d. Calculate induced V(t) along the cable section and total voltage (sum).
    #         e. DSTL voltage calculation
    #         f. Total voltage calculation

    #         Preprocessings:
    #         ---------------
    #         1. Read mapping file
    #         2. To Ben's csv file to BEZpy reaable .txt
    #         """
    #         logger.info(f"Load cable section mappings")
    #         logger.info(f"Start simulation....")
    #         synb = SynB(self.args)
    #         synb.run()
    #         cab = Cable(self.args, None, synb.Bfield, synb.components)
    #         cab.run_nodal_analysis()
    # #         e = EventAnalysis(o, self.out_dirs["base"], self.verbose,
    # #                           self.syne, self.E_syn, self.out_dirs["base"])
    # #         e.calclulate_total_parameters()
    #         logger.info(f"Simulation end!")
    #         return

    @staticmethod
    def run(args):
        """
        This is a static method that helps to run
        the simulation step-by-step manner
        """
        sim = Simulation(args)
        if args.opcode == 0:
            sim.B_field_simulation()
        elif args.opcode == 1:
            sim.analyze_event_data()
        return


# Script run can also be done via main program
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-sch", "--schema_folder", default="input/schemas/", help="Schema folder"
    )
    parser.add_argument(
        "-prm",
        "--param_file",
        default="input/jsons/SynB.json",
        help="Simulation Parameter File",
    )
    parser.add_argument(
        "-opc", "--opcode", default=0, help="Operation code [0/1]", type=int
    )
    args = parser.parse_args()
    logger.info("Parameter list for simulation ")
    utils.print_rec(args)
    isValid, o = utils.validate_jsons(args.param_file, args.schema_folder, args.opcode)
    if isValid:
        logger.info(f"Fork folders to save data/figs/params.")
        utils.create_structures(o, args.param_file)
        utils.print_rec(o)
        Simulation.run(o)
    logger.info(f"Simulation end! Clear local files.")
    os.system("rm -rf `find -type d -name .ipynb_checkpoints`:")
    os.system("rm -rf __pycache__/")
