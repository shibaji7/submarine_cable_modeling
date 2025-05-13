"""uqanalysis.py: Module is used to implement all simulations for UQ"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import numpy as np


class OceanEarthProfile(object):
    def __init__(self):
        self.segments = np.arange(9).astype(str)
        # Create
        return


if __name__ == "__main__":
    OceanEarthProfile()
