import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scubas.sources.current_sources import Current


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


def plot_real_imag_plots(c, d):
    setups()

    fig, axes = plt.subplots(nrows=2, ncols=2, dpi=240, figsize=(10, 6))
    axes = axes.ravel()
    ax = axes[0]
    ax.plot(
        d / 1e3,
        1e9 * np.real(c.B_fields[0]["total"]),
        ls="-",
        color="k",
        lw=1.2,
        label="Total",
    )
    ax.plot(
        d / 1e3,
        1e9 * np.real(c.B_fields[0]["external"]),
        ls="--",
        color="r",
        lw=0.8,
        label="External",
    )
    ax.plot(
        d / 1e3,
        1e9 * np.real(c.B_fields[0]["internal"]),
        ls="--",
        color="b",
        lw=0.8,
        label="Internal",
    )
    ax.set_xlim(d[0] / 1e3, d[-1] / 1e3)
    ax.legend(loc=1)

    ax = axes[1]
    ax.plot(d / 1e3, 1e9 * np.real(c.B_fields[1]["total"]), ls="-", color="k", lw=1.2)
    ax.plot(
        d / 1e3, 1e9 * np.real(c.B_fields[1]["external"]), ls="--", color="r", lw=0.8
    )
    ax.plot(
        d / 1e3, 1e9 * np.real(c.B_fields[1]["internal"]), ls="--", color="b", lw=0.8
    )
    ax.set_xlim(d[0] / 1e3, d[-1] / 1e3)

    ax = axes[2]
    ax.plot(d / 1e3, 1e9 * np.imag(c.B_fields[0]["total"]), ls="-", color="k", lw=1.2)
    ax.plot(
        d / 1e3, 1e9 * np.imag(c.B_fields[0]["external"]), ls="--", color="r", lw=0.8
    )
    ax.plot(
        d / 1e3, 1e9 * np.imag(c.B_fields[0]["internal"]), ls="--", color="b", lw=0.8
    )
    ax.set_xlim(d[0] / 1e3, d[-1] / 1e3)

    ax = axes[3]
    ax.plot(d / 1e3, 1e9 * np.imag(c.B_fields[1]["total"]), ls="-", color="k", lw=1.2)
    ax.plot(
        d / 1e3, 1e9 * np.imag(c.B_fields[1]["external"]), ls="--", color="r", lw=0.8
    )
    ax.plot(
        d / 1e3, 1e9 * np.imag(c.B_fields[1]["internal"]), ls="--", color="b", lw=0.8
    )
    ax.set_xlim(d[0] / 1e3, d[-1] / 1e3)

    axes[0].set_ylabel(r"$B_{real}$, nT")
    axes[2].set_ylabel(r"$B_{imag}$, nT")
    axes[2].set_xlabel("Distance, km")
    axes[3].set_xlabel("Distance, km")

    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("1998_B.png", bbox_inches="tight")

    fig, axes = plt.subplots(nrows=2, ncols=1, dpi=240, figsize=(5, 6))
    axes = axes.ravel()

    ax = axes[0]
    ax.plot(d / 1e3, 1e3 * np.real(c.E_fields[0]["total"]), ls="-", color="r", lw=0.8)
    ax.set_xlim(d[0] / 1e3, d[-1] / 1e3)
    ax.set_ylim(-6, 1)
    ax = axes[1]
    ax.plot(d / 1e3, 1e3 * np.imag(c.E_fields[0]["total"]), ls="-", color="r", lw=0.8)
    ax.set_xlim(d[0] / 1e3, d[-1] / 1e3)
    ax.set_ylim(-6, 1)

    axes[0].set_ylabel(r"$E_{real}$, V/km")
    axes[1].set_ylabel(r"$E_{imag}$, V/km")
    axes[1].set_xlabel("Distance, km")

    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("1998_E.png", bbox_inches="tight")
    return


def test_current_sources():
    T = np.array([5]) * 60
    omega = 2 * np.pi / T
    d = np.linspace(-1000, 1000, 2001)
    c = Current(h=100, a=200, I=1, omega=omega, layer=0, return_Z=True)
    # print(c.Z, c.Z_out, c.p)
    # c.p = c.site.calcZ(1/T, layer=0)[1]/(1j*omega)
    # print(c.p, c.omega, omega, c.freq, 1/T)
    c.compute_B_field(d)

    plot_real_imag_plots(c, d)
    return


if __name__ == "__main__":
    test_current_sources()
