import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import numpy as np
import pandas as pd
from scipy import constants as C

caseB = False

if caseB:
    # s2 >> s1
    sigma1, sigma2 = 1/3000, 1/3
    t1 = 100000
    flim, M = [1e-8, 1e0], 10
    freqs = np.linspace(flim[0], flim[1], int(flim[1] / (M*flim[0])) + 1)
    fig = plt.figure(dpi=300, figsize=(2,4))
    ax = fig.add_subplot(211)
    r = np.sqrt(sigma1/sigma2)
    g = np.sqrt(1j*2*np.pi*freqs*C.mu_0*sigma1)*t1
    Q = (r+np.tanh(g))/(1+r*np.tanh(g))
    ax.semilogx(freqs, np.abs(np.tanh(g)), "r", lw=1.0, ls="-", label="tanh")
    ax.set_ylabel(r"$|\tanh{\frac{\tau_1}{\delta_A}}|$", color="r")
    ax.text(0.7, 0.9, r"Case B [$\sigma1<<\sigma2$]"+"\n"+fr"$\tau_1=100.0$ km", ha="right", va="center", transform=ax.transAxes)
    ax = ax.twinx()
    ax.semilogx(freqs, np.angle(np.tanh(g), deg=True), "b", lw=1.0, ls="-")
    ax.set_yticks([0, 30, 45, 60, 90])
    ax.set_ylabel(r"$\angle \tanh{\frac{\tau_1}{\delta_A}}$ $[^\circ]$", color="b")
    ax.set_ylim(0, 90)
    ax.set_xticks([1e-8, 1e-4, 1e0])
    ax.set_xticks([])
    ax.set_xlim(flim)


    ax = fig.add_subplot(212)
    ax.semilogx(freqs, np.abs(Q), "r", lw=1.0, ls="-", label="tanh")
    ax.set_ylabel(r"$|Q|$", color="r")
    ax = ax.twinx()
    ax.semilogx(freqs, np.angle(Q, deg=True), "b", lw=1.0, ls="-")
    ax.set_yticks([0, 30, 45, 60, 90])
    ax.set_ylabel(r"$\angle Q$ $[^\circ]$", color="b")
    ax.set_ylim(0, 90)
    ax.set_xticks([1e-8, 1e-4, 1e0])
    ax.set_xlim(flim)

    fig.subplots_adjust(wspace=0.05,hspace=0.05)
    fig.savefig("cmnt4.png", bbox_inches="tight")
else:
    # s1 >> s2
    flim, M = [1e-8, 1e0], 5
    freqs = np.linspace(flim[0], flim[1], int(flim[1] / (M*flim[0])) + 1)

    sigma1, sigma2 = 1/3, 1/3000
    t1 = 10000

    fig = plt.figure(dpi=300, figsize=(2,4))
    ax = fig.add_subplot(211)
    r = np.sqrt(sigma1/sigma2)
    g = np.sqrt(1j*2*np.pi*freqs*C.mu_0*sigma1)*t1
    Q = (r+np.tanh(g))/(1+r*np.tanh(g))
    ax.semilogx(freqs, np.abs(1./np.tanh(g)), "r", lw=1.0, ls="-", label="tanh")
    ax.set_ylabel(r"$|\coth{\frac{\tau_1}{\delta_A}}|$", color="r")
    ax.text(0.7, 0.9, r"Case A [$\sigma1>>\sigma2$]"+"\n"+fr"$\tau_1=10.0$ km", ha="right", va="center", transform=ax.transAxes)
    ax = ax.twinx()
    ax.semilogx(freqs, np.angle(1./np.tanh(g), deg=True), "b", lw=1.0, ls="-")
    ax.set_yticks([-90, -60, -45, -30, 0])
    ax.set_ylabel(r"$\angle \coth{\frac{\tau_1}{\delta_A}}$ $[^\circ]$", color="b")
    ax.set_ylim(-90, 0)
    ax.set_xticks([1e-8, 1e-4, 1e0])
    ax.set_xticks([])
    ax.set_xlim(flim)


    ax = fig.add_subplot(212)
    ax.semilogx(freqs, np.abs(Q), "r", lw=1.0, ls="-", label="tanh")
    ax.set_ylabel(r"$|Q|$", color="r")
    ax = ax.twinx()
    ax.semilogx(freqs, np.angle(Q, deg=True), "b", lw=1.0, ls="-")
    ax.set_yticks([-90, -60, -45, -30, 0])
    ax.set_ylabel(r"$\angle Q$ $[^\circ]$", color="b")
    ax.set_ylim(-90, 0)
    ax.set_xticks([1e-8, 1e-4, 1e0])
    ax.set_xlim(flim)

    fig.subplots_adjust(wspace=0.05,hspace=0.05)
    fig.savefig("cmnt5.png", bbox_inches="tight")