import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
import matplotlib.pyplot as plt

G      = 6.674e-11
k      = 1.380649e-23
Msolar = 1.989e30
mH     = 1.6726e-27
c      = 3e8
GeV_to_kg = 1e9 * 1.602e-19 / c**2
pc  = 3.0857e16
kpc = 3.0857e19

class NFW_profile:
    def __init__(self, rho0, rs):
        self.rho0, self.rs = rho0, rs

    def density(self, r):
        x = r / self.rs
        return self.rho0 / (x * (1 + x)**2)

    def velocity_disp(self, d, n=2000):
        r = np.linspace(1e-6, d, n)
        M = np.trapezoid(4 * np.pi * r**2 * self.density(r), r)
        return np.sqrt(G * M / d)   # v_c(d)

    def DM_mass_local(self, R, d):
        return self.density(d) * (4/3) * np.pi * R**3

def jeans_energy(Rdm, Mnm, Rnm, halo, d, T):
    Mdm = halo.DM_mass_local(Rdm, d)
    vsq = halo.velocity_disp(d) ** 2

    NM_self = -3 * G * Mnm**2 / (5 * Rnm)
    DM_self = -3 * G * Mdm**2 / (5 * Rdm)
    int_energy = (-1.5 * G * Mnm * Mdm / Rdm
                  + 0.3 * G * Mnm * Mdm * Rnm**2 / Rdm**3)
    KE = 0.5 * Mdm * vsq + 1.5 * (Mnm/mH) * k * T

    return NM_self + DM_self + int_energy + KE


def bound_fraction(v_esc, sigma):
    pdf = lambda v: np.sqrt(2/np.pi) * v**2 / sigma**3 * np.exp(-v**2/(2*sigma**2))
    upper = min(v_esc, 50*sigma)
    frac, _ = quad(pdf, 0, upper)
    return 1.0 if v_esc > 50*sigma else frac


halo = NFW_profile(rho0=0.51 * GeV_to_kg / (1e-2)**3, rs=8.1 * kpc)

Mnm = 10 * Msolar     # baryon clump mass
Rnm = 3.05e16           # cm
T   = 150               # gas temperature [K]

dist = np.logspace(-1, 5) * pc
bound_dm = []
valid_dist = []
R_max = []
for d in dist:
    try:
        Rmax = brentq(jeans_energy, Rnm * 1.01, 1e5 * pc, args=(Mnm, Rnm, halo, d, T))
        Mdm_naive = halo.DM_mass_local(Rmax, d)
        sigma = halo.velocity_disp(d) / np.sqrt(3)
        v_esc = np.sqrt(2 * G * (Mnm) / Rmax)
        frac = bound_fraction(v_esc, sigma)
        bound_dm.append(Mdm_naive * frac)
        valid_dist.append(d)
        R_max.append(Rmax)
    except:
        continue
    
valid_dist = np.array(valid_dist)
bound_dm = np.array(bound_dm)
R_max = np.array(R_max)

fig, ax = plt.subplots(1, 2)

ax[0].loglog(valid_dist / pc, bound_dm / Msolar)
ax[0].set_xlabel("Distance from Halo Centre (pc)")
ax[0].set_ylabel(r"Bounded DM ($M_\odot$)")

ax[1].loglog(valid_dist/pc, R_max/pc)
ax[1].set_xlabel("Distance from Halo Centre (pc)")
ax[1].set_ylabel("Collapsed DM radius")

plt.show()
"""
plt.figure(figsize=(4, 3))
plt.loglog(valid_dist/pc, bound_dm/Msolar)
plt.xlabel("Distance from Halo Centre (pc)")
plt.ylabel(r"Bounded DM ($M_\odot$)")
plt.tight_layout()
plt.savefig("./img/accretedDM.pgf")
plt.show()
"""