import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def setups(size=12):
    plt.style.use(["science", "ieee"])
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [
        "Tahoma",
        "DejaVu Sans",
        "Lucida Grande",
        "Verdana",
    ]
    mpl.rcParams.update(
        {"xtick.labelsize": size, "ytick.labelsize": size, "font.size": size}
    )
    return


def plot_real_imag_plots(dist, B, file_name=""):
    setups()
    fig, axes = plt.subplots(nrows=1, ncols=1, dpi=240, figsize=(5, 3))
    ax = axes
    ax.plot(dist / 1e3, 1e9 * np.real(B), ls="-", color="k", lw=1.2, label="Real")
    ax.plot(dist / 1e3, 1e9 * np.imag(B), ls="-", color="r", lw=0.8, label="Imag")
    ax.set_xlim(dist[-1] / 1e3, dist[0] / 1e3)
    ax.legend(loc=1)
    ax.set_ylabel("B, nT")
    ax.set_xlabel("Distance, km")
    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig(file_name, bbox_inches="tight")
    return
