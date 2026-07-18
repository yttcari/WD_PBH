from WhiteDwarf.wd_setup import WhiteDwarf
import WhiteDwarf.source as source
import numpy as np
from tqdm import tqdm
import pandas as pd
from WhiteDwarf.constants import *
import os

def iterate(rhoc_scaled, source_function, max_iter=20, epsilon=5e-4, DEBUG=False, **kwargs):
    wd = WhiteDwarf(Ye=0.5, rhoc_scaled=rhoc_scaled, source=source_function, Z=6, **kwargs)
    wd.hydro_integrate(DEBUG=DEBUG)
    para0 = wd.thermo_integrate(DEBUG=DEBUG)
    print(f"r (km) {wd.rbar2r(para0[0]):.3e}, rho (g/cc): {wd.rhobar2rho(para0[1]):.3e}, M (Msolar): {wd.mbar2m(para0[2]):.3e}, T (K): {(para0[3] * wd.rho0 * c ** 2/ a) ** 0.25:.3e}, Luminosity: {source_function.luminosity(r=wd.R_profile[0] * wd.R0, m=wd.M_profile[-1] * wd.M0):.3e}")

    print("Finished initial setup.\n")
    for n in range(max_iter):
        print(f"Round {n}: hydro")
        wd.hydro_integrate(DEBUG=DEBUG)
        print(f"Round {n}: thermal")
        para = wd.thermo_integrate(DEBUG=DEBUG)
        print(f"r (km) {wd.rbar2r(para[0]):.3e}, rho (g/cc): {wd.rhobar2rho(para[1]):.3e}, M (Msolar): {wd.mbar2m(para[2]):.3e}, T (K): {(para[3] * wd.rho0/ a) ** 0.25:.3e}, Luminosity: {source_function.luminosity(r=wd.R_profile[0] * wd.R0, m=wd.M_profile[-1] * wd.M0):.3e}")
        

        eps = np.nanmax(np.abs(para-para0)/para0)
        if (eps <= epsilon).all():
            print(f"Epsilon: {np.mean(eps):.3e}, Converges, BREAK")
            return wd
            
        else:
            print(f"Finished iter: {n}, epsilon: {np.mean(eps):.3e} CONTINUE")
            para0 = para

    print("===== DOESN'T CONVERGE =====")

PBH_list = [1e15]
rho_list = [1e-1, 1e0, 1e1, 1e2, 1e3]
N_list = np.logspace(5, 10, 10)
for N in N_list:
    for PBH_M in PBH_list:
        PBH = source.singlePBH(PBH_M, N=N)
        new_dir = os.path.join(PBH.source_type, f'PBH_M{PBH_M:.1e}_N{N:.1e}/')
        os.makedirs(new_dir, exist_ok=True)
        for rho in tqdm(rho_list):
            print(f"Working on rho: {rho}, M: {PBH_M}")
            output_path = os.path.join(new_dir, f'rho0{rho:.3e}.csv')
            if os.path.exists(output_path):
                print(f"Skipping rho={rho:.3e}, file already exists")
                continue
            wd = iterate(rho, PBH, DEBUG=False, dr=1e-3, epsilon=1e-4, max_iter=20)
            if wd is not None:
                print(f"White dwarf with PBH: {wd.rbar2r(wd.R_profile[-1]):.3e} km, {wd.mbar2m(wd.M_profile[-1]):.3e} Msolar")
                wd.write_csv(output_path)

"""
for PBH_M in PBH_list:
    PBH = source.uniformPBH(PBH_M, N=1e10, R=100*km2cm)
    new_dir = os.path.join(PBH.source_type, f'PBH_M{PBH_M:.1e}/')
    os.makedirs(new_dir, exist_ok=True)
    for rho in tqdm(rho_list):
        print(f"Working on rho: {rho}, M: {PBH_M}")
        output_path = os.path.join(new_dir, f'rho0{rho:.3e}.csv')
        if os.path.exists(output_path):
            print(f"Skipping rho={rho:.3e}, file already exists")
            continue
        wd = iterate(rho, PBH, DEBUG=False, dr=1e-3, epsilon=1e-4, max_iter=20)
        if wd is not None:
            print(f"White dwarf with PBH: {wd.rbar2r(wd.R_profile[-1]):.3e} km, {wd.mbar2m(wd.M_profile[-1]):.3e} Msolar")
            wd.write_csv(output_path)

for PBH_M in PBH_list:
    PBH = source.gaussianPBH(PBH_M, N0=10, sigma=100*km2cm)
    new_dir = os.path.join(PBH.source_type, f'PBH_M{PBH_M:.1e}/')
    os.makedirs(new_dir, exist_ok=True)
    for rho in tqdm(rho_list):
        print(f"Working on rho: {rho}, M: {PBH_M}")
        output_path = os.path.join(new_dir, f'rho0{rho:.3e}.csv')
        if os.path.exists(output_path):
            print(f"Skipping rho={rho:.3e}, file already exists")
            continue
        wd = iterate(rho, PBH, DEBUG=False, dr=1e-3, epsilon=1e-4, max_iter=20)
        if wd is not None:
            print(f"White dwarf with PBH: {wd.rbar2r(wd.R_profile[-1]):.3e} km, {wd.mbar2m(wd.M_profile[-1]):.3e} Msolar")
            wd.write_csv(output_path)
"""