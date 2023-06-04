# Import modules
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import re
import scipy

# Source CRYODIR
sys.path.append(os.environ["CRYODIR"])
# User defined functions
from postProcessing import concat_of_fobjs


def source_p_exp(
    filename,
    exp_path=os.environ["CRYODIR"] + "/validation/non_isobaric/Seo_2010_LN2/exp_data/",
):
    df_exp = pd.read_csv(exp_path + filename)
    t_exp = df_exp["Time_min"].values * 60
    p_exp = df_exp["Pressure_Pa"].values
    return t_exp, p_exp


import subprocess

# Source exp data
t_exp, p_exp = source_p_exp("Seo_LF70_2.5W_exp.csv")

# Source model data


def eta_u_fobj(x):
    """
    Objective function to minimize the error in
    pressure, vapour temperature and liquid temperature

    x: np.array
        x[0]: evaporative fraction eta_w
        x[1]: overall heat transfer coefficient U
    """
    # Get current working directory
    cwd = os.getcwd()

    # Read transport properties to get the substring
    # to replace the value of eta_evap

    U_names = ["U_L", "U_dry", "U_roof", "U_bot"]

    # Create empty list of U values
    subs_U = []
    U_lines = []
    pattern = re.compile("eta_evap")

    # Store contents of the file
    with open(cwd + "/constant/transportProperties", "r") as file:
        f = file.read()

    with open(cwd + "/constant/transportProperties", "r") as file:
        # Store transport properties file
        for line in file:
            for match in re.finditer(pattern, line):
                eta_line = line
                eta = line.split()[1]

            for U_name in U_names:
                for match in re.finditer(U_name, line):
                    U_lines.append(line)
                    subs_U.append(line.split()[1])

    # Replace new guess of eta_evap
    # debug prints
    # print(eta_line)
    # print(subs_U)
    # print(U_lines)

    new_eta_line = eta_line.replace(eta, str(x[0]) + ";")

    new_U_lines = []
    # Update overall heat transfer coefficients
    for idx, U_line in enumerate(U_lines):
        new_U_lines.append(U_line.replace(subs_U[idx], str(x[1]) + ";"))

    print(new_eta_line)
    print(new_U_lines)
    # print(subs_eta)

    # Writing to file
    with open(cwd + "/constant/transportProperties", "w") as file:
        # Replace eta_evap
        f = f.replace(eta_line, new_eta_line)
        for idx, U_line in enumerate(U_lines):
            f = f.replace(U_line, new_U_lines[idx])

        file.write(f)

        # Create liquid temperature files
    with open("del.log", "w") as outfile:
        subprocess.run(["rm", "-r", "postProcessing"], stdout=outfile)

    # Run OpenFOAM!

    with open("simulation.log", "w") as outfile:
        subprocess.run("buoyantBoussinesqPvapFoam", stdout=outfile)

    # Create liquid temperature files
    with open("tl.log", "w") as outfile:
        subprocess.run(["postProcess", "-func", "singleGraph"], stdout=outfile)

    # Read pressure function object
    df = concat_of_fobjs(cwd, fobj_name="p_average")

    # Extract time and pressure
    t_mod = df["Time"].values
    p_mod = df["p_average"].values

    p_mod_int = np.interp(t_exp, t_mod, p_mod)
    RMSE_p = np.mean(np.sqrt((p_mod_int - p_exp) ** 2))
    print("Pressure RMSE = %.1f Pa" % RMSE_p)

    # Get temperature
    # Plot and rescale z on the basis of CFD mesh
    height = 0.213  # m
    LF = 0.717

    exp_path = os.environ["CRYODIR"] + "/validation/non_isobaric/Seo_2010_LN2/exp_data/"

    # Source experimental data
    df = pd.read_csv(exp_path + "LF70_2.5W_Tz.csv")

    T_exp = []
    for i in range(3, 6):
        T_exp.append(df.loc[df["phase"] == "vapour"].iloc[:, i].values)

    tl_exp = []
    for i in range(3, 6):
        tl_exp.append(df.loc[df["phase"] == "liquid"].iloc[:, i].values)

    zl_exp = df.loc[df["phase"] == "liquid"].iloc[:, 1].values / 1000

    z_exp = df.loc[df["phase"] == "vapour"].iloc[:, 1].values / 1000
    # Translate the profile to the itnerface
    z_exp = z_exp - height * LF

    temp_df = pd.read_csv(cwd + "/2400/vap_temp.csv", skiprows=2, nrows=100)
    z_df = pd.read_csv(cwd + "/2400/vap_mesh.csv", skiprows=2, nrows=100)

    z = z_df.iloc[:, 0].values
    temp_mod = temp_df.iloc[:, 0].values

    RMSE_tv = np.sqrt(np.mean((np.interp(z_exp, z, temp_mod) - T_exp[1]) ** 2))
    print("Vapour temperature RMSE = %.3f K " % RMSE_tv)

    tl_cfd_name = cwd + "/postProcessing/singleGraph/2400/line_T.xy"
    tl_cfd_half = pd.read_csv(tl_cfd_name, sep="\t", header=None)
    zl_cfd = tl_cfd_half.iloc[:, 0].values
    tl_cfd = tl_cfd_half.iloc[:, 1].values

    RMSE_tl = np.sqrt(np.mean((np.interp(zl_exp, zl_cfd, tl_cfd) - tl_exp[1]) ** 2))
    print("Liquid temperature RMSE = %.3f K " % RMSE_tl)

    L = (
        70 * RMSE_p / np.mean(p_exp)
        + 20 * RMSE_tv / np.mean(T_exp[1])
        + 10 * RMSE_tl / np.mean(tl_exp[1])
    )
    print("Loss = %.3f, eta = %.3f, U = %.3f" % (L, x[0], x[1]))

    return L


# U, eta
eta_0 = 0.97
U_0 = 0.087
x0 = [eta_0, U_0]

# Lower and upper bounds for parameter fitting
lb = [0.95, 0.07]
ub = [1, 0.15]

scipy.optimize.least_squares(eta_u_fobj, x0, bounds=(lb, ub), max_nfev=14)

# eta_u_fobj([0.96, 0.087])
