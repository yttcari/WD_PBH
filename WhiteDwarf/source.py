import numpy as np
from WhiteDwarf.constants import *
import os
import pandas as pd

hbar = h / (2 * np.pi)

class proton_decay:
    def __init__(self, k, **kwargs):
        self.k = k # s-1
    
    def luminosity(self, r, m, **kwargs):
        """
        m: enclosed mass in g
        k: decay constant in 1/s
        E_gamma: gamma ray energy in erg
        mp: proton mass in g
        """
        L = m * self.k * proton_energy/ mp # unit: erg/s
        return L

class No_source:
    def __init__(self, **kwargs):
        pass

    def get_mass(self):
        return 0
    
    def luminosity(self, **kwargs):
        return 0

class singlePBH:
    def __init__(self, pbhM, N=1, **kwargs):
        self.pbhM = pbhM
        self.N = N
        print(f"Hawking temperature of the PBH is {hbar * c ** 3 / (8 * np.pi * G * pbhM * kB):.3e}")
        self.source_type = 'SinglePBH'

    def get_mass(self):
        # return mass for TOV calculation
        return self.pbhM * self.N

    def hawking_T(self):
        return hbar * c ** 3 / (8 * np.pi * G * self.pbhM * kB)

    def get_spectra(self, spectra_type):
        mass = f"{self.pbhM:.0e}"
        spectra_dir = f"/Users/caritsang/Library/CloudStorage/OneDrive-TheChineseUniversityofHongKong/01 Project/FYP_Sem1/blackhawk_v2.3/results/M{mass}/instantaneous_{spectra_type}_spectra.txt"

        if not os.path.exists:
            print("The entered PBH mass has no seconary spectra computed. Please compute it using BlackHawk before you run this. TERMINATE")
            raise ValueError("No Secondary Spectra found")
        else:
            df = pd.read_csv(spectra_dir, sep='\s+', skiprows=1)
            return df

    def luminosity(self, **kwargs):
        #return hbar * c**6 / (15360 * np.pi * G**2 * self.pbhM**2) * self.N
        GEV_TO_ERG = 1.602176634e-3      # 1 GeV = 1.602e-3 erg

        sec_df = self.get_spectra('secondary')
        pri_df = self.get_spectra('primary')

        # primary spectra
        def add_pri():
            dN_dE = pri_df['photon'].values / GEV_TO_ERG
            E = pri_df['energy/particle'].values * GEV_TO_ERG
            lumin = np.trapz(dN_dE, E)

            return lumin
    
        L = add_pri()

        # secondary spectra
        E = sec_df['energy/particle'].values * GEV_TO_ERG

        sec_species_list = ["photon", "electron", "proton", "neutron"]
        for sp in sec_species_list:
            # Convert dN/dE from GeV^-1 → erg^-1
            dN_dE = sec_df[sp].values / GEV_TO_ERG

            # Integrate over energy
            lumin = np.trapz(dN_dE, E)
            L += lumin

        return L * self.N

class uniformPBH:
    def __init__(self, pbhM, N, R, **kwargs):
        self.pbhM = pbhM # mass of PBH [g]
        self.N = N # total number of PBH
        self.R = R # truncate the PBH ball [cm]
        self.source_type = 'UniformPBH'

    def get_mass(self, r):
        # r in cm
        return self.pbhM * self.N * min(r / self.R, 1) ** 3
    
    def get_spectra(self, spectra_type):
        mass = f"{self.pbhM:.0e}"
        spectra_dir = f"/Users/caritsang/Library/CloudStorage/OneDrive-TheChineseUniversityofHongKong/01 Project/FYP_Sem1/blackhawk_v2.3/results/M{mass}/instantaneous_{spectra_type}_spectra.txt"

        if not os.path.exists:
            print("The entered PBH mass has no seconary spectra computed. Please compute it using BlackHawk before you run this. TERMINATE")
            raise ValueError("No Secondary Spectra found")
        else:
            df = pd.read_csv(spectra_dir, sep='\s+', skiprows=1)
            return df

    def luminosity(self, r, **kwargs):
        #return hbar * c**6 / (15360 * np.pi * G**2 * self.pbhM**2) * self.N
        GEV_TO_ERG = 1.602176634e-3      # 1 GeV = 1.602e-3 erg

        sec_df = self.get_spectra('secondary')
        pri_df = self.get_spectra('primary')

        # primary spectra
        def add_pri():
            dN_dE = pri_df['photon'].values / GEV_TO_ERG
            E = pri_df['energy/particle'].values * GEV_TO_ERG
            lumin = np.trapz(dN_dE, E)

            return lumin
    
        L = add_pri()

        # secondary spectra
        E = sec_df['energy/particle'].values * GEV_TO_ERG

        sec_species_list = ["photon", "electron", "proton", "neutron"]
        for sp in sec_species_list:
            # Convert dN/dE from GeV^-1 → erg^-1
            dN_dE = sec_df[sp].values / GEV_TO_ERG

            # Integrate over energy
            lumin = np.trapz(dN_dE, E)
            L += lumin

        return L * (self.N * min(r / self.R, 1) ** 3)

    
class gaussianPBH:
    def __init__(self, pbhM, N0, sigma, **kwargs):
        self.pbhM  = pbhM    # mass of individual PBH [g]
        self.N0    = N0      # central number
        self.sigma = sigma   # spread [cm]
        self.source_type = 'GaussianPBH'

    def _N_enclosed(self, r):
        from scipy.special import erf
        r   = np.asarray(r)
        u   = r / (np.sqrt(2) * self.sigma)
        return 4 * np.pi * self.N0 * self.sigma**3 * (
            np.sqrt(np.pi / 2) * erf(u) - (r / self.sigma) * np.exp(-u**2)
        )

    def get_mass(self, r):
        return self.pbhM * self._N_enclosed(r)

    def get_spectra(self, spectra_type):
        mass = f"{self.pbhM:.0e}"
        spectra_dir = f"/Users/caritsang/Library/CloudStorage/OneDrive-TheChineseUniversityofHongKong/01 Project/FYP_Sem1/blackhawk_v2.3/results/M{mass}/instantaneous_{spectra_type}_spectra.txt"

        if not os.path.exists:
            print("The entered PBH mass has no seconary spectra computed. Please compute it using BlackHawk before you run this. TERMINATE")
            raise ValueError("No Secondary Spectra found")
        else:
            df = pd.read_csv(spectra_dir, sep='\s+', skiprows=1)
            return df

    def luminosity(self, r, **kwargs):
        #return hbar * c**6 / (15360 * np.pi * G**2 * self.pbhM**2) * self.N
        GEV_TO_ERG = 1.602176634e-3      # 1 GeV = 1.602e-3 erg

        sec_df = self.get_spectra('secondary')
        pri_df = self.get_spectra('primary')

        # primary spectra
        def add_pri():
            dN_dE = pri_df['photon'].values / GEV_TO_ERG
            E = pri_df['energy/particle'].values * GEV_TO_ERG
            lumin = np.trapz(dN_dE, E)

            return lumin
    
        L = add_pri()

        # secondary spectra
        E = sec_df['energy/particle'].values * GEV_TO_ERG

        sec_species_list = ["photon", "electron", "proton", "neutron"]
        for sp in sec_species_list:
            # Convert dN/dE from GeV^-1 → erg^-1
            dN_dE = sec_df[sp].values / GEV_TO_ERG

            # Integrate over energy
            lumin = np.trapz(dN_dE, E)
            L += lumin

        return L * self._N_enclosed(r)