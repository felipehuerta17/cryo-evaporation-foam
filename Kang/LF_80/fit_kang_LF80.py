# Parameter optimisation routine
# This needs to be executed from the case directory

# Import modules
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy
import pickle

# Source CRYODIR
sys.path.append(os.environ["CRYODIR"])

# User defined functions
from postProcessing import concat_of_fobjs
from optimise import source_p_exp, source_p_mod
from optimise.data_source import source_tl_exp, source_tv_exp
from optimise.fit_non_isobaric import eta_u_fobj

# Define experimental path and filename
exp_path = os.environ["CRYODIR"] + "/validation/non_isobaric/Kang_2018_LN2/exp_data/"

# WORKSPACE VARIABLES
# Initial liuid filling
LF = 0.83
# Tank height
height = 0.758054  #

# Source experimental pressure
t_exp, p_exp = source_p_exp("LF_80_P.csv", exp_path)

# Experimental temperature timesteps
times = ["1200", "2400", "3600", "4800", "5220"]
# Debug
# times = ["1"]

# Source experimental liquid temperature
zl_exp, tl_exp = source_tl_exp("LF_80_Tvz.csv", exp_path, height, LF)

# Source experimental vapour temperature
zv_exp, tv_exp = source_tv_exp("LF_80_Tvz.csv", exp_path, height, LF)

# Source full path of the simulation case
casedir = (
    os.environ["CRYODIR"]
    + "OpenFOAM/non_isobaric/singlephase/Kang/neq_new_etawall/LF_80_opt3_U028"
)

# U, eta
eta_0 = 0.98
U_0 = 0.1
x0 = [eta_0, U_0]
# Lower and upper bounds for parameter fitting

# Evaporative coefficient
eta_range = (0.98, 0.984)

# Overall heat transfer coefficient
U_range = (0.39, 0.43)

# Number of mesh points to evaluate
n_mesh = 3

# Brute force global optimisation
x0, fval, grid, Jout = scipy.optimize.brute(
    eta_u_fobj,
    args=(t_exp, p_exp, zl_exp, tl_exp, zv_exp, tv_exp, height, LF, times, casedir),
    ranges=(eta_range, U_range),
    Ns=n_mesh,
    finish=None,
    full_output=True,
)

# Save into dictionary
opti_dict = {"x0": x0, "fval": fval, "grid": grid, "Jout": Jout}

with open("optimal.pkl", "wb") as f:
    pickle.dump(opti_dict, f)

print(opti_dict)

# Dual annealing
# res = scipy.optimize.dual_annealing(
#     eta_u_fobj, x0=x0, bounds=((0.96, 1), (0.06, 0.10)), maxiter=6, maxfun=20,
# )
# print(res)

# eta_u_fobj([0.96, 0.087])
