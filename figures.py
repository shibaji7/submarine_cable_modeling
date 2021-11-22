"""figures.py: Module is used to implement plot figures for publications"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
os.system("rm -rf `find -type d -name .ipynb_checkpoints`:")
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])

from models import OceanModel

def create_stack_plots_TFs():
    bins = ["data/OceanModels/RhoZ_Bin%d.txt"%x for x in np.arange(1,10)]
    tfs = []
    for f in np.arange(1,10):
        tf = OceanModel.getOceanModel(f).get_TFs()
        tfs.append(tf)
        print(" Bins %02d"%f)
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=180, figsize=(3,3))
    colors = ["darkred", "darkblue", "darkgreen"]*3
    linestyles = ["-"]*3 + ["--"]*3 + ["-."]*3
    for f in np.arange(1,10):
        ax.loglog(tfs[f-1].freq, np.absolute(tfs[f-1].Ef2Bs), colors[f-1], ls=linestyles[f-1],
                  lw=0.8, label=r"$\text{Bin}_{%02d}$"%f)
    ax.legend(loc=0, prop={"size": 8}, bbox_to_anchor=(1.01, 1.))
    ax.set_xlabel(r"$f_0$, (Hz)")
    ax.set_ylabel(r"$|\frac{E_f}{B_s}|, (mV/km/nT)$")
    ax.set_ylim([1e-2,1e0])
    ax.set_xlim(1e-4,1e-2)
    fig.savefig("docs/TF.stack.png", bbox_inches="tight")
    return

if __name__ == "__main__":
    create_stack_plots_TFs()