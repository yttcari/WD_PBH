"""
Writes a BlackHawk parameters.txt file from the variables below, then
runs ./BlackHawk_tot.x on it. Run this script from the BlackHawk root
directory (the one containing BlackHawk_tot.x).

$ python write_param.py
"""

import numpy as np
import subprocess
import sys
from pathlib import Path
import os
import shutil

params = {
	"destination_folder": "PBH_evolution/M1e+15",
	"full_output": 0,
	"interpolation_method": 0,

	"metric": 0,

	"BH_number": 1,
	"Mmin": 1e+15,
	"Mmax": 1e+15,
	"param_number": 1,
	"amin": 0.,
	"amax": 0.,
	"Qmin": 0.,
	"Qmax": 0.,

	"epsilon_LQG": 1.5,
	"a0_LQG": 0.,
	"n": 0.,

	"spectrum_choice": 0,
	"spectrum_choice_param": 0,

	"amplitude_lognormal": 1.,
	"amplitude_lognormal2": 1.,
	"stand_dev_lognormal": 1.,
	"crit_mass_lognormal": 1.,

	"amplitude_powerlaw": 1.,
	"eqstate_powerlaw": 0.3333,

	"amplitude_critical_collapse": 1.,
	"crit_mass_critical_collapse": 1.,

	"amplitude_uniform": 1.,

	"stand_dev_param_Gaussian": 1.,
	"mean_param_Gaussian": 0.5,

	"table": "spin_distribution_BH.txt",

	"tmin_manual": 0,
	"tmin": 1.e-30,
	"limit": 5000,
	"BH_remnant": 0,
	"M_remnant": 1e-4,

	"E_number": 10,
	"Emin": 1e-6,
	"Emax": 5,

	"grav": 1,
	"add_DM": 0,
	"m_DM": 1.,
	"spin_DM": 0.,
	"dof_DM": 1.,

	"primary_only": 0,

	"hadronization_choice": 3,

	"bondi_accretion": 0,
	"cs": 2e8,
	"lambda_bondi": 1,
	"rho_bg": 1e6,

	"unruh_accretion": 1,
	"T_unruh": 1e6,
	"m_unruh": 938.272e-3 * 6,  # evaluated in Python: the C parser reads a single number, not an expression
	"rho_unruh": 1e6,
}

def sweep_param(WD_rho, MPBH, BONDI, UNRUH, params):
	params['destination_folder'] = f"PBH_evolution/bondi/M{MPBH:.2e}/rho{WD_rho:.2e}"
	params['Mmax'] = MPBH
	params['Mmin'] = MPBH

	params['bondi_accretion'] = BONDI
	params['unruh_accretion'] = UNRUH
	
	params['rho_bg'] = WD_rho
	params['rho_unruh'] = WD_rho

	return params

def write_parameters_file(path, params):
	lines = []
	for key, value in params.items():
		lines.append(f"{key} = {value}")
	Path(path).write_text("\n".join(lines) + "\n")

def run_blackhawk(param_file, executable="./BlackHawk_tot.x"):
	#subprocess.run(["make", "clean"])
	#subprocess.run(["make", executable[2:-2]])
	if not Path(executable).exists():
		sys.exit(f"ERROR: '{executable}' not found. Run this script from the BlackHawk root directory.")
	subprocess.run([executable, param_file])

if __name__ == "__main__":
    os.chdir('../../blackhawk_v2.3')
    
    param_file = "./param.txt"
    rho_list = np.logspace(6, 9, 15)
    MPBH_list = np.logspace(13, 20, 15)
    
    for rho in rho_list:
        for MPBH in MPBH_list:
            params = sweep_param(WD_rho=rho, MPBH=MPBH, BONDI=1, UNRUH=0, params=params)
            os.makedirs(f"./results/{params['destination_folder']}", exist_ok=True)
            actual_output_dir = os.path.join("results", params['destination_folder'])
            if os.path.exists(actual_output_dir):
                shutil.rmtree(actual_output_dir)
                
            write_parameters_file(path=param_file, params=params)
            print(f"Running BlackHawk for M={MPBH:.3e}, rho={rho:.3e}")
            run_blackhawk(param_file)