import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import numpy as np
import pandas as pd
from scipy import constants as C

def calc_Z(f, rho1, rho2, l1):
    sigma1, sigma2 = 1/rho1, 1/rho2
    i = 1j
    k1 = np.sqrt(i*2*np.pi*f*C.mu_0*sigma1)
    eta1 = i*2*np.pi*f/k1
    K2 = np.sqrt(i*2*np.pi*f/(C.mu_0*sigma2))
    alpha = 2*k1*l1
    header = (K2 * (1+np.exp(-alpha))) + (eta1 * (1-np.exp(-alpha)))
    footer = (K2 * (1-np.exp(-alpha))) + (eta1 * (1+np.exp(-alpha)))
    Eq31 = eta1 * (header/footer) /(1.0e-3 / C.mu_0)
    uni_rho2 = K2 /(1.0e-3 / C.mu_0)
    uni_rho1 = np.sqrt(i*2*np.pi*f/(C.mu_0*sigma1)) /(1.0e-3 / C.mu_0)
    return dict(Eq31=Eq31, l1=l1, rho1=rho1, rho2=rho2, uni_rho2=uni_rho2,uni_rho1=uni_rho1)

def abline():
    gca = plt.gca()
    gca.set_autoscale_on(False)
    gca.plot(gca.get_xlim(),gca.get_ylim(), ls="-.", lw=0.4, color="k")
    return

def abxline():
    gca = plt.gca()
    gca.set_autoscale_on(False)
    gca.plot(gca.get_xlim(),[4e-1,4e-1], ls="-.", lw=0.4, color="k")
    return

flim, M = [1e-6, 1e0], 1
f = np.linspace(flim[0], flim[1], int(flim[1] / (M*flim[0])) + 1)
rho1, rho2 = 3, 3000
l1 = 10000
case1Z = calc_Z(f, rho1, rho2, l1)
case2Z = calc_Z(f, rho2, rho1, l1*10)


fig = plt.figure(dpi=300, figsize=(2,4))
ax = fig.add_subplot(211)
ax.loglog(f, np.abs(case1Z["Eq31"]), "r", lw=1.0, ls="-", label="non-U")
ax.loglog(f, np.abs(case1Z["uni_rho1"]), "r", lw=0.6, ls="--", zorder=3, label=r"U($\rho$=3 $\Omega.m$)")
ax.loglog(f, np.abs(case1Z["uni_rho2"]), "r", lw=0.6, ls=":", zorder=3, label=r"U($\rho$=3000 $\Omega.m$)")
ax.legend(loc='upper left', bbox_to_anchor=(1.1, 1.0),)
ax.set_ylabel(r"Amplitude [mV/km/nT]", color="r")
ax.set_xlim(flim)
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_ylim(1e-3,1e3)
abxline()
ax.set_yticks([1e-3, 1e0, 1e3])
ax.text(0.9, 0.9, fr"Case A, $\tau_0={l1/1000}$ km", ha="right", va="center", transform=ax.transAxes)
ax = ax.twinx()
ax.semilogx(f, np.angle(case1Z["Eq31"], deg=True), "b", lw=1.0, ls="-")
ax.semilogx(f, np.angle(case1Z["uni_rho1"], deg=True), "b", lw=0.6, ls="--", zorder=3)
ax.semilogx(f, np.angle(case1Z["uni_rho2"], deg=True), "b", lw=0.6, ls=":", zorder=3)
ax.set_yticks([0, 30, 45, 60, 90])
ax.set_ylabel(r"Phase $[^\circ]$", color="b")
ax.set_ylim(0, 90)
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_xlim(flim)


ax = fig.add_subplot(212)
ax.loglog(f, np.abs(case2Z["Eq31"]), "r", lw=1.0, ls="-", label="non-U")
ax.loglog(f, np.abs(case2Z["uni_rho2"]), "r", lw=0.6, ls="--", zorder=3, label=r"U($\rho$=3 $\Omega.m$)")
ax.loglog(f, np.abs(case2Z["uni_rho1"]), "r", lw=0.6, ls=":", zorder=3, label=r"U($\rho$=3000 $\Omega.m$)")
ax.set_xlabel(r"Frequency [Hz]")
ax.set_ylabel(r"Amplitude [mV/km/nT]", color="r")
ax.set_xlim(flim)
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_ylim(1e-3,1e3)
ax.set_yticks([1e-3, 1e0, 1e3])
abline()
ax.text(0.9, 0.1, fr"Case B, $\tau_0={l1*10/1000}$ km", ha="right", va="center", transform=ax.transAxes)
ax = ax.twinx()
ax.semilogx(f, np.angle(case2Z["Eq31"], deg=True), "b", lw=1.0, ls="-")
ax.semilogx(f, np.angle(case2Z["uni_rho2"], deg=True), "b", lw=0.6, ls="--", zorder=3)
ax.semilogx(f, np.angle(case2Z["uni_rho1"], deg=True), "b", lw=0.6, ls=":", zorder=3)
ax.set_yticks([0, 30, 45, 60, 90])
ax.set_ylabel(r"Phase $[^\circ]$", color="b")
ax.set_ylim(0, 90)
ax.set_xlim(flim)
ax.set_xticks([1e-6, 1e-3, 1e0])


fig.savefig("cmnt0.png", bbox_inches="tight")



flim, M = [1e-6, 1e0], 1
f = np.linspace(flim[0], flim[1], int(flim[1] / (M*flim[0])) + 1)
rho1, rho2 = 3, 3000
l1 = 10000
k1 = np.sqrt(1j*2*np.pi*f*C.mu_0/rho1)
K2 = np.sqrt(1j*2*np.pi*f/(C.mu_0/rho2))
eta1 = 1j*2*np.pi*f/k1
B = ( (K2*(1+np.exp(-2*k1*l1))) + (eta1*(1-np.exp(-2*k1*l1))) )/((eta1*(1+np.exp(-2*k1*l1))) + (K2*(1-np.exp(-2*k1*l1))))
index = 1000

fig = plt.figure(dpi=300, figsize=(4,4))
ax = fig.add_subplot(221)
ax.loglog(f, np.abs(2*k1*l1), "r", lw=1.0, ls="-")
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_ylim([1e-3, 1e6])
ax.set_ylabel(r"Amplitude [mV/km/nT]", color="r")
ax.set_xlim(flim)
ax = ax.twinx()
ax.semilogx(f, np.angle(2*k1*l1, deg=True), "b", lw=1.0, ls="-")
ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax.set_ylim(-90, 90)
ax.set_xlim(flim)
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_yticklabels([])
ax.text(0.9, 0.9, r"(a) $2k_1 l_1$", ha="right", va="center", transform=ax.transAxes)
ax.text(0.1, 1.1, r"$\rho_1=3$ $\Omega.m$, $\rho_2=3000$ $\Omega.m$" + "\n" + r"$l_1=10.0$ km",ha="left", va="center", transform=ax.transAxes)
ax.axvline(f[index], color="k", ls="--", lw=1.2)
ax.text(0.1, 0.1, r"$\angle\frac{1}{k_1}$(1 mHz)=%.2f"%np.angle(1./k1, deg=True)[index],
         ha="left", va="center", transform=ax.transAxes)
r1 = (eta1 - K2) /(eta1 + K2)
print(2*k1[index]*l1, np.abs(2*k1[index]*l1), np.angle(2*k1[index]*l1, deg=True))
print(-np.exp(-2*k1[index]*l1), np.abs(-np.exp(-2*k1[index]*l1)), 
      np.angle(-np.exp(-2*k1[index]*l1), deg=True))
print(eta1[index], np.abs(eta1[index]), 
      np.angle(eta1[index], deg=True))
print(K2[index], np.abs(K2[index]), 
      np.angle(K2[index], deg=True))
print(r1[index])
print(1-r1[index]*np.exp(-2*k1[index]*l1), np.abs(1-r1[index]*np.exp(-2*k1[index]*l1)),
      np.angle(1-r1[index]*np.exp(-2*k1[index]*l1), deg=True))
print(1+r1[index]*np.exp(-2*k1[index]*l1), np.abs(1+r1[index]*np.exp(-2*k1[index]*l1)),
      np.angle(1+r1[index]*np.exp(-2*k1[index]*l1), deg=True))
print(1/(1+r1[index]*np.exp(-2*k1[index]*l1)), np.abs(1/(1+r1[index]*np.exp(-2*k1[index]*l1))),
      np.angle(1/(1+r1[index]*np.exp(-2*k1[index]*l1)), deg=True))
print((1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1)), np.abs((1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1))),
      np.angle((1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1)), deg=True))
print(eta1[index]*(1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1)),
      np.abs(eta1[index]*(1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1))),
      np.angle(eta1[index]*(1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1)), deg=True))

ax = fig.add_subplot(222)
ax.loglog(f, np.abs(eta1), "r", lw=1.0, ls="-")
ax.set_yticklabels([])
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_ylim([1e-3, 1e6])
ax.set_xlim(flim)
ax = ax.twinx()
ax.semilogx(f, np.angle(eta1, deg=True), "b", lw=1.0, ls="-")
ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax.set_ylabel(r"Phase $[^\circ]$", color="b")
ax.set_ylim(-90, 90)
ax.set_xlim(flim)
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.text(0.9, 0.9, r"(b) $\eta_1$", ha="right", va="center", transform=ax.transAxes)
ax.axvline(f[index], color="k", ls="--", lw=1.2)
ax.text(0.9, 1.05, "$f_0$=1 mHz", ha="right", va="center", transform=ax.transAxes)
ax.text(0.1, 0.1, r"$\angle\eta_1$(1 mHz)=%.2f"%np.angle(eta1, deg=True)[index],
         ha="left", va="center", transform=ax.transAxes)

ax = fig.add_subplot(223)
ax.loglog(f, np.abs(K2), "r", lw=1.0, ls="-")
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_ylim([1e-3, 1e6])
ax.set_ylabel(r"Amplitude [mV/km/nT]", color="r")
ax.set_xlabel(r"Frequency [Hz]")
ax.set_xlim(flim)
ax = ax.twinx()
ax.semilogx(f, np.angle(K2, deg=True), "b", lw=1.0, ls="-")
ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax.set_ylim(-90, 90)
ax.set_xlim(flim)
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_yticklabels([])
ax.text(0.9, 0.9, r"(c) $K_2$", ha="right", va="center", transform=ax.transAxes)
ax.axvline(f[index], color="k", ls="--", lw=1.2)
ax.text(0.1, 0.1, r"$\angle K_2$(1 mHz)=%.2f"%np.angle(K2, deg=True)[index],
         ha="left", va="center", transform=ax.transAxes)

ax = fig.add_subplot(224)
ax.loglog(f, np.abs(B), "r", lw=1.0, ls="-")
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_ylim([1e-3, 1e6])
ax.set_yticklabels([])
ax.set_xlabel(r"Frequency [Hz]")
ax.set_xlim(flim)
ax = ax.twinx()
ax.semilogx(f, np.angle(B, deg=True), "b", lw=1.0, ls="-")
ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax.set_ylabel(r"Phase $[^\circ]$", color="b")
ax.set_ylim(-90, 90)
ax.set_xlim(flim)
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.text(0.9, 0.9, r"(d) $B$", ha="right", va="center", transform=ax.transAxes)
ax.text(0.1, 0.1, r"$\angle B$(1 mHz)=%.2f"%np.angle(B, deg=True)[index],
         ha="left", va="center", transform=ax.transAxes)
ax.axvline(f[index], color="k", ls="--", lw=1.2)
fig.savefig("cmnt0.png", bbox_inches="tight")


rho1, rho2 = 3000, 3
l1 = 100000
k1 = np.sqrt(1j*2*np.pi*f*C.mu_0/rho1)
K2 = np.sqrt(1j*2*np.pi*f/(C.mu_0/rho2))
eta1 = 1j*2*np.pi*f/k1
B = ( (K2*(1+np.exp(-2*k1*l1))) + (eta1*(1-np.exp(-2*k1*l1))) )/((eta1*(1+np.exp(-2*k1*l1))) + (K2*(1-np.exp(-2*k1*l1))))
r1 = (eta1 - K2) /(eta1 + K2)
print("\n\n")
print(2*k1[index]*l1, np.abs(2*k1[index]*l1), np.angle(2*k1[index]*l1, deg=True))
print(-np.exp(-2*k1[index]*l1), np.abs(-np.exp(-2*k1[index]*l1)), 
      np.angle(-np.exp(-2*k1[index]*l1), deg=True))
print(eta1[index], np.abs(eta1[index]), 
      np.angle(eta1[index], deg=True))
print(K2[index], np.abs(K2[index]), 
      np.angle(K2[index], deg=True))
print(">>", r1[index])
print(1-r1[index]*np.exp(-2*k1[index]*l1), np.abs(1-r1[index]*np.exp(-2*k1[index]*l1)),
      np.angle(1-r1[index]*np.exp(-2*k1[index]*l1), deg=True))
print(1+r1[index]*np.exp(-2*k1[index]*l1), np.abs(1+r1[index]*np.exp(-2*k1[index]*l1)),
      np.angle(1+r1[index]*np.exp(-2*k1[index]*l1), deg=True))
print(1/(1+r1[index]*np.exp(-2*k1[index]*l1)), np.abs(1/(1+r1[index]*np.exp(-2*k1[index]*l1))),
      np.angle(1/(1+r1[index]*np.exp(-2*k1[index]*l1)), deg=True))
print((1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1)), np.abs((1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1))),
      np.angle((1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1)), deg=True))
print(eta1[index]*(1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1)),
      np.abs(eta1[index]*(1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1))),
      np.angle(eta1[index]*(1-r1[index]*np.exp(-2*k1[index]*l1))/(1+r1[index]*np.exp(-2*k1[index]*l1)), deg=True))


fig = plt.figure(dpi=300, figsize=(4,4))
ax = fig.add_subplot(221)
ax.loglog(f, np.abs(1./k1), "r", lw=1.0, ls="-")
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_ylim([1e-3, 1e9])
ax.set_ylabel(r"Amplitude [mV/km/nT]", color="r")
ax.set_xlim(flim)
ax = ax.twinx()
ax.semilogx(f, np.angle(1./k1, deg=True), "b", lw=1.0, ls="-")
ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax.set_ylim(-90, 90)
ax.set_xlim(flim)
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_yticklabels([])
ax.axvline(f[index], color="k", ls="--", lw=1.2)
ax.text(0.9, 0.9, r"(a) $k_1$", ha="right", va="center", transform=ax.transAxes)
ax.text(0.1, 1.1, r"$\rho_1=3000$ $\Omega.m$, $\rho_2=3$ $\Omega.m$" + "\n" + r"$l_1=100.0$ km",ha="left", va="center", transform=ax.transAxes)
ax.text(0.1, 0.1, r"$\angle\frac{1}{k_1}$(1 mHz)=%.2f"%np.angle(1./k1, deg=True)[index],
         ha="left", va="center", transform=ax.transAxes)

ax = fig.add_subplot(222)
ax.loglog(f, np.abs(eta1), "r", lw=1.0, ls="-")
ax.set_yticklabels([])
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_ylim([1e-3, 1e9])
ax.set_xlim(flim)
ax = ax.twinx()
ax.semilogx(f, np.angle(eta1, deg=True), "b", lw=1.0, ls="-")
ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax.set_ylabel(r"Phase $[^\circ]$", color="b")
ax.set_ylim(-90, 90)
ax.set_xlim(flim)
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.axvline(f[index], color="k", ls="--", lw=1.2)
ax.text(0.9, 0.9, r"(b) $\eta_1$", ha="right", va="center", transform=ax.transAxes)
ax.text(0.9, 1.05, "$f_0$=1 mHz", ha="right", va="center", transform=ax.transAxes)
ax.text(0.1, 0.1, r"$\angle\eta_1$(1 mHz)=%.2f"%np.angle(eta1, deg=True)[index],
         ha="left", va="center", transform=ax.transAxes)

ax = fig.add_subplot(223)
ax.loglog(f, np.abs(K2), "r", lw=1.0, ls="-")
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_ylim([1e-3, 1e9])
ax.set_ylabel(r"Amplitude [mV/km/nT]", color="r")
ax.set_xlabel(r"Frequency [Hz]")
ax.set_xlim(flim)
ax = ax.twinx()
ax.semilogx(f, np.angle(K2, deg=True), "b", lw=1.0, ls="-")
ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax.set_ylim(-90, 90)
ax.set_xlim(flim)
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_yticklabels([])
ax.axvline(f[index], color="k", ls="--", lw=1.2)
ax.text(0.9, 0.9, r"(c) $K_2$", ha="right", va="center", transform=ax.transAxes)
ax.text(0.1, 0.1, r"$\angle K_2$(1 mHz)=%.2f"%np.angle(K2, deg=True)[index],
         ha="left", va="center", transform=ax.transAxes)

ax = fig.add_subplot(224)
ax.loglog(f, np.abs(B), "r", lw=1.0, ls="-")
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.set_ylim([1e-3, 1e9])
ax.set_yticklabels([])
ax.set_xlabel(r"Frequency [Hz]")
ax.set_xlim(flim)
ax = ax.twinx()
ax.semilogx(f, np.angle(B, deg=True), "b", lw=1.0, ls="-")
ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax.set_ylabel(r"Phase $[^\circ]$", color="b")
ax.set_ylim(-90, 90)
ax.set_xlim(flim)
ax.set_xticks([1e-6, 1e-3, 1e0])
ax.text(0.9, 0.9, r"(d) $B$", ha="right", va="center", transform=ax.transAxes)
ax.axvline(f[index], color="k", ls="--", lw=1.2)
ax.text(0.1, 0.1, r"$\angle B$(1 mHz)=%.2f"%np.angle(B, deg=True)[index],
         ha="left", va="center", transform=ax.transAxes)
fig.savefig("cmnt0.png", bbox_inches="tight")