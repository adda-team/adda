# -*- coding: utf-8 -*-
"""
V5.0

@author: Mattia Andrini, mattia.andrini@unicatt.it

Computes integrated Forward/Backward scattering efficiencies over hemispherical domains, 
accounting for full azimuthal (phi) dependence.

Features:
- Computes integrated Forward/Backward scattering efficiencies (Csca) over user-defined hemispherical domains
- User-friendly GUI for simulation setup.
- Automated ADDA simulations over a wavelength spectrum.
- Fully adjustable Theta-Phi angular grid.
- Native support for substrate mode.
- Support for custom particle shapes and OpenCL (GPU) acceleration.
"""

import time
import numpy as np
import pylab
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import shutil
from pathlib import Path

# --- 1. SETUP PATHS & GLOBAL VARIABLES ---

# Determine the directory where this script is located
basedir = Path(__file__).resolve().parent
Mainfolder = basedir
SimResfolder = basedir / "SimulationResults"

# Ensure the results folder exists immediately
SimResfolder.mkdir(parents=True, exist_ok=True)

# Global variables for data storage
wavelength = None
PS_n = None
PS_k = None
sub_n = None
ADDA_EXE_PATH = None
ADDA_OCL_EXE_PATH = None

# --- 2. ENVIRONMENT & EXECUTABLE SETUP ---

def setup_adda_environment():
    """
    Finds adda.exe, adda_ocl.exe, and required DLLs using relative paths logic.
    Prioritizes the local folder, then searches up the directory tree.
    Copies missing files to the script folder to ensure runnability.
    """
    global ADDA_EXE_PATH, ADDA_OCL_EXE_PATH
    
    # List of files required for ADDA to run (Standard + OpenCL version)
    required_files = ["adda.exe", "adda_ocl.exe", "clFFT.dll", "libfftw3-3.dll"]
    
    # Potential paths where files might be, relative to this script.
    search_dirs = [
        basedir,
        basedir / "win64",
        basedir.parent / "win64",
        basedir.parent.parent / "src" / "win64",
        basedir.parent.parent / "win64"
    ]

    all_found = True
    
    for filename in required_files:
        file_found = False
        local_file = basedir / filename
        
        # Check potential directories for this specific file
        for folder in search_dirs:
            potential_path = folder / filename
            
            if potential_path.exists():
                file_found = True
                
                # Update global paths for executables
                if filename == "adda.exe":
                    ADDA_EXE_PATH = potential_path.resolve()
                    print(f"[INFO] Found ADDA executable at: {ADDA_EXE_PATH}")
                elif filename == "adda_ocl.exe":
                    ADDA_OCL_EXE_PATH = potential_path.resolve()
                    print(f"[INFO] Found ADDA OCL executable at: {ADDA_OCL_EXE_PATH}")
                
                # Copy to local folder if found elsewhere
                if potential_path != local_file:
                    try:
                        if not local_file.exists():
                            print(f"[INFO] Copying {filename} to script folder...")
                            shutil.copy(potential_path, local_file)
                            
                            # Update executable path to local if we just copied it
                            if filename == "adda.exe":
                                ADDA_EXE_PATH = local_file
                            elif filename == "adda_ocl.exe":
                                ADDA_OCL_EXE_PATH = local_file
                    except Exception as e:
                        print(f"[WARN] Could not copy {filename} locally: {e}")
                
                break # Stop searching for this file once found
        
        if not file_found:
            # Warning about adda_ocl.exe absence, don't fail the whole setup if it's missing,
            # as standard adda.exe is the critical one for CPU runs.
            if filename == "adda_ocl.exe":
                print(f"[WARNING] Optional file 'adda_ocl.exe' not found. GPU modes will not work.")
            else:
                print(f"[WARNING] Could not find required file: {filename}")
                all_found = False

    if not all_found:
        messagebox.showwarning("Missing ADDA Files", 
                               "Could not find 'adda.exe' or required DLLs (clFFT.dll, libfftw3-3.dll) in standard locations.\n\n"
                               "Please manually place them in the same folder as this script.")
        return False
        
    return True

# --- 3. HANDLE scat_params.dat file ---
#sca_params.dat is created or overwritten accounting for the selected GUI angular resolution before each simulation.

def create_scat_params(save_directory, theta_n, theta_min, theta_max, phi_n, phi_min, phi_max):
    filepath = save_directory / "scat_params.dat"
    file_content = f"""# Description of a set of angles for which to calculate Mueller matrix
global_type=grid
N=2
pairs=
0 0
30 90

theta:
type=range
N={theta_n}
min={theta_min}
max={theta_max}
values=
0
10
20
30

phi:
type=range
N={phi_n}
min={phi_min}
max={phi_max}
values=
0
90

phi_integr:
min=0
max=360
Jmin=2
Jmax=5
eps=0
equiv=true
periodic=true
"""
    try:
        with open(filepath, 'w') as f:
            f.write(file_content)
    except IOError as e:
        print(f"Error writing to file {filepath}: {e}")

def generate_bat_file(wavelength, PS_n, PS_k, sub_n, k, a): #builds the batch file sent to ADDA.
    bat_filename = f"run_{k}.bat"
    
    wavelength_micro = float(wavelength) / 1000
    lambda_value = int(wavelength)
    eq_rad_value = float(a) / 1000
    dist_surf_value = float(distance_from_surf) / 1000
    
    run_dir_name = f"run{lambda_value}_{particle_shape}_sub{substrate_presence}_g{grid_value}_m{PS_n}_R{str(a)}_{Flag}"

    with open(bat_filename, 'w') as f:
        f.write("@echo off\n")
        
        if GPU_acc == "no":
            command = f'"{ADDA_EXE_PATH}" -scat_matr muel -lambda {wavelength_micro} -m {PS_n} {PS_k} '
        else:
            command = f'"{ADDA_OCL_EXE_PATH}" -gpu {GPU_acc} -scat_matr muel -lambda {wavelength_micro} -m {PS_n} {PS_k} '

        if size_type == "Equivalent radius":
            command += f"-eq_rad {eq_rad_value} "
        elif size_type == "Size Along X axis":
            command += f"-size {eq_rad_value} "

        if substrate_presence == "yes":
            command += f"-surf {dist_surf_value} {sub_n} 0 -prop 0 0 -1 "
        else:
            command += "-prop 0 0 1 "

        if custom_shape == "no":
            command += f"-shape {particle_shape} -grid {grid_value} "
        else:
            command += f"-shape read {particle_shape}.geom "

        if Csca_option == "yes":
            command += "-Csca "
        
        if store_scat_grid_option == "yes":
            command += "-store_scat_grid "
        else:
            command += f"-ntheta {ntheta_value} "

        command += f"-dir {run_dir_name} \n"
        f.write(command)

    return bat_filename

def run_batch_file(bat_file, k):
    start_time = time.time()
    try:
        subprocess.run(bat_file, shell=True, check=True)
        print(f"{wavelength[k]} simulation took: {time.time() - start_time:.2f} seconds")
    except subprocess.CalledProcessError as e:
        print(f"Error running batch file: {e}")

def go(a): #starts ADDA simulation
    if not ADDA_EXE_PATH:
        if not setup_adda_environment(): return

    if store_scat_grid_option == "yes":
        create_scat_params(Mainfolder, theta_n=(ntheta_value + 1), theta_min=0, 
                          theta_max=str(int(theta_stop)), phi_n=(len(phi_grid)), 
                          phi_min=0, phi_max=str(int(phi_stop)))

    start_time = time.time()
    bat_files = []
    
    for k in range(wavelength.shape[0]):
        bat_file = generate_bat_file(wavelength[k], PS_n[k], PS_k[k], sub_n[k], k, a)
        bat_files.append((bat_file, k))

    print("Starting parallel simulations...")
    with ThreadPoolExecutor(max_workers=parallel_processes) as executor:
        futures = [executor.submit(run_batch_file, bat_file, k) for bat_file, k in bat_files]
        for future in as_completed(futures):
            future.result()
            
    if keep_bat_files == "no":
        for bat_file, _ in bat_files:
            try:
                os.remove(bat_file)
            except Exception as e:
                print(f"Could not delete {bat_file}: {e}")
                
    print(f"R={a}nm simulation took: {time.time() - start_time:.2f} seconds")

def kill_folder(a):
    for j, wl in enumerate(wavelength):
        run = f"{int(wl):03d}"
        folder_name = f"run{run}_{particle_shape}_sub{substrate_presence}_g{number_of_dipoles_in_a_row}_m{str(PS_n[j])}_R{str(a)}_{Flag}"
        folder = Mainfolder / folder_name
        if folder.exists():
            shutil.rmtree(folder)
            print(f'Folder deleted: {folder}')
        else:
            print('Next folder...')

# --- 4. DATA PROCESSING & MIE FUNCTIONS ---

def QSca_Mie_for(m, lambda0, a):
    theta = np.arange(0, 0.5*np.pi+1e-8, np.pi/100000)
    x = 2*np.pi*a/lambda0
    try:
        import miepython
        S1, S2 = miepython._mie_S1_S2(m, x, np.cos(theta), 6)
        return x**(-2) * np.sum((np.abs(S1)**2 + np.abs(S2)**2) * np.sin(theta) * np.pi/100000)
    except: return 0.0

def QSca_Mie_back(m, lambda0, a):
    theta = np.arange(0.5*np.pi, np.pi+1e-8, np.pi/100000)
    x = 2*np.pi*a/lambda0
    try:
        import miepython
        S1, S2 = miepython._mie_S1_S2(m, x, np.cos(theta), 6)
        return x**(-2) * np.sum((np.abs(S1)**2 + np.abs(S2)**2) * np.sin(theta) * np.pi/100000)
    except: return 0.0

def Qext_Mie(m, lambda0, a):
    x = 2*np.pi*a/lambda0
    try:
        import miepython
        S1, S2 = miepython._mie_S1_S2(m, x, np.array([1.0, 1.01]), 6)
        return 4 * (1/x**2) * S1[0].real
    except: return 0.0

# --- 5. UPLOAD & PLOT FUNCTIONS ---

def upload(a):  # Upload values from a simulation, calculates scattering quantities and saves them in a .txt file when the phi scattering option is disabled.
    # Handle the run parameters
    if wavelength.shape[0] > 1:
        initial_run = int(wavelength[0])
        step_lambda = int(wavelength[1] - wavelength[0])
        final_run = initial_run + step_lambda * wavelength.shape[0]
    else:
        initial_run = int(wavelength[0])
        step_lambda = 1  # no step if only one point
        final_run = initial_run + step_lambda  # function discontinued

    S11_muel = np.zeros((len(theta_grid), wavelength.shape[0]))  # Just initialize the arrays

    if Csca_option == "yes":
        Pol_Y = np.zeros((6, wavelength.shape[0]))
        Pol_X = np.zeros((6, wavelength.shape[0]))
    if Csca_option == "no":
        Pol_Y = np.zeros((4, wavelength.shape[0]))
        Pol_X = np.zeros((4, wavelength.shape[0]))

    Q_ext_ADDA_autoY = np.zeros(wavelength.shape[0])
    Q_abs_ADDA_autoY = np.zeros(wavelength.shape[0])
    Q_ext_ADDA_autoX = np.zeros(wavelength.shape[0])
    Q_abs_ADDA_autoX = np.zeros(wavelength.shape[0])
    Q_sca_ADDA_intY = np.zeros(wavelength.shape[0])
    Q_sca_ADDA_intX = np.zeros(wavelength.shape[0])
    C_ext_ADDA_autoY = np.zeros(wavelength.shape[0])
    C_abs_ADDA_autoY = np.zeros(wavelength.shape[0])
    C_ext_ADDA_autoX = np.zeros(wavelength.shape[0])
    C_abs_ADDA_autoX = np.zeros(wavelength.shape[0])
    C_sca_ADDA_intY = np.zeros(wavelength.shape[0])
    C_sca_ADDA_intX = np.zeros(wavelength.shape[0])

    for j, wl in enumerate(wavelength):
        i = int(wl)  # wavelength value, used for folder naming

        if i < 10:
            run = "00" + str(i)
        if 10 <= i < 100:
            run = "0" + str(i)
        else:
            run = str(i)

        filename = f"run{run}_{particle_shape}_sub{substrate_presence}_g{number_of_dipoles_in_a_row}_m{str(PS_n[j])}_R{str(a)}_{Flag}"

        MuelMatrx = Mainfolder / filename / "mueller"  # location of mueller matrix

        # j = int((i - initial_run)/step_lambda)

        if wavelength.shape[0] == 1:
            theta = np.loadtxt(MuelMatrx, delimiter=" ", skiprows=1, usecols=0)[:len(theta_grid)]

        if i == initial_run:
            theta = np.loadtxt(MuelMatrx, delimiter=" ", skiprows=1, usecols=0)[:len(theta_grid)]
        S11_muel[:len(theta), j] = np.loadtxt(MuelMatrx, delimiter=" ", skiprows=1, usecols=1)[:len(theta)]

        CrossSecY = Mainfolder / filename / "CrossSec-Y"  # Y cross sections

        CrossSecX_candidate = Mainfolder / filename / "CrossSec-X"  # X cross sections

        # If CrossSec-X exists, use it. Otherwise, fall back to CrossSec-Y (CrossSec-Y may be the only output due to symmetry)
        if os.path.exists(CrossSecX_candidate):
            CrossSecX = CrossSecX_candidate
        else:
            CrossSecX = CrossSecY

        columnY = []
        with open(CrossSecY) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 3:
                    columnY.append([float(parts[2])])

        columnX = []
        with open(CrossSecX) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 3:
                    columnX.append([float(parts[2])])

        Pol_Y[:, j] = np.array(columnY).flatten()
        Pol_X[:, j] = np.array(columnX).flatten()

        # Loads cross sections and efficiency factors from the CrossSec files, for each wavelength index j.

        Q_abs_ADDA_autoY[j] = Pol_Y[3, j]
        Q_ext_ADDA_autoY[j] = Pol_Y[1, j]
        Q_abs_ADDA_autoX[j] = Pol_X[3, j]
        Q_ext_ADDA_autoX[j] = Pol_X[1, j]

        C_abs_ADDA_autoY[j] = Pol_Y[2, j]
        C_ext_ADDA_autoY[j] = Pol_Y[0, j]
        C_abs_ADDA_autoX[j] = Pol_X[2, j]
        C_ext_ADDA_autoX[j] = Pol_X[0, j]

        if Csca_option == "yes":
            Q_sca_ADDA_intY[j] = Pol_Y[5, j]  # these output are shown only if -Csca option is selected and it is useful to check that everything is correct.
            Q_sca_ADDA_intX[j] = Pol_X[5, j]

            C_sca_ADDA_intY[j] = Pol_Y[4, j]
            C_sca_ADDA_intX[j] = Pol_X[4, j]
        else:
            Q_sca_ADDA_intY[j] = 0
            Q_sca_ADDA_intX[j] = 0

            C_sca_ADDA_intY[j] = 0
            C_sca_ADDA_intX[j] = 0

    Q_abs_ADDA_auto = (Q_abs_ADDA_autoY + Q_abs_ADDA_autoX) / 2  # Averages for different directions
    Q_ext_ADDA_auto = (Q_ext_ADDA_autoY + Q_ext_ADDA_autoX) / 2
    Q_sca_ADDA_int = (Q_sca_ADDA_intY + Q_sca_ADDA_intX) / 2

    C_abs_ADDA_auto = (C_abs_ADDA_autoY + C_abs_ADDA_autoX) / 2
    C_ext_ADDA_auto = (C_ext_ADDA_autoY + C_ext_ADDA_autoX) / 2
    C_sca_ADDA_int = (C_sca_ADDA_intY + C_sca_ADDA_intX) / 2

    theta_90_lab = int(np.round(len(theta_grid) / 2) + 1)
    theta_91_lab_lab = int(np.round(len(theta_grid) / 2) + 2)

    print(f'90° is set as row {theta_90_lab}')  # prints the row of the mueller file assigned to theta=90°. Useful to check consistency.

    if substrate_presence == "no":
        theta_angle_for = theta[:theta_90_lab]  # Defines set of angles for forward scattering (0-90°).
        theta_angle_back = theta[theta_91_lab_lab:]
    if substrate_presence == "yes":  # while working in substrate mode, theta angle must be reversed to assume that positive z axis direction is the impinging light direction.
        theta_angle_for = theta[theta_91_lab_lab:]
        theta_angle_back = theta[:theta_90_lab]

    Q_ext_ADDA = np.zeros(len(wavelength))
    Q_sca_for_ADDA = np.zeros(len(wavelength))
    Q_sca_back_ADDA = np.zeros(len(wavelength))
    C_ext_ADDA = np.zeros(len(wavelength))
    C_sca_for_ADDA = np.zeros(len(wavelength))
    C_sca_back_ADDA = np.zeros(len(wavelength))
    Q_sca_Mie_for = np.zeros(len(wavelength))
    Q_sca_Mie_back = np.zeros(len(wavelength))
    Q_ext_Mie = np.zeros(len(wavelength))
    C_sca_ADDA_int_tot = np.zeros(len(wavelength))

    x = np.zeros(len(wavelength))

    for k in range(len(wavelength)):
        x[k] = 2 * np.pi * a / wavelength[k]  # Compute size parameter

        if substrate_presence == "no":
            S11_muel_for = S11_muel[:theta_90_lab, k]
            S11_muel_back = S11_muel[theta_91_lab_lab:, k]

        if substrate_presence == "yes":  # while working in substrate mode, theta angle must be reversed to assume that positive z axis direction is the impinging light direction.
            S11_muel_for = S11_muel[theta_91_lab_lab:, k]
            S11_muel_back = S11_muel[:theta_90_lab:, k]

        if substrate_presence == "no":  # in substrate mode, the scattering efficiency of the light propagating in the substrate must be divided by the refractive index of the substrate squared.
            sub_n[k] = 1

        dtheta = theta_grid[1] - theta_grid[0]  # Takes the theta step from the mueller file.

        if size_type == "Size Along X axis":
            Q_sca_for_ADDA[k] = (np.pi * (a / 1000) ** 2 / particle_cross_section) * 4 * 2 * x[k] ** (-2) * np.sum(
                S11_muel_for * np.sin(theta_angle_for * np.pi / 180) * dtheta * 2 * np.pi / 360) / sub_n[k] ** 2  # Integrates forward scattering efficiencies assuming SPHERES!
            Q_sca_back_ADDA[k] = (np.pi * (a / 1000) ** 2 / particle_cross_section) * 4 * 2 * x[k] ** (-2) * np.sum(
                S11_muel_back * np.sin(theta_angle_back * np.pi / 180) * dtheta * 2 * np.pi / 360)
            Q_ext_ADDA[k] = (np.pi * (a / 1000) ** 2 / particle_cross_section) * 4 * (4 / x[k] ** 2) * (S11_muel_for[0])  # Compute Qext
        else:
            Q_sca_for_ADDA[k] = (np.pi * (a / 1000) ** 2 / particle_cross_section) * 2 * x[k] ** (-2) * np.sum(
                S11_muel_for * np.sin(theta_angle_for * np.pi / 180) * dtheta * 2 * np.pi / 360) / sub_n[k] ** 2  # Integrates forward scattering efficiencies assuming SPHERES!
            Q_sca_back_ADDA[k] = (np.pi * (a / 1000) ** 2 / particle_cross_section) * 2 * x[k] ** (-2) * np.sum(
                S11_muel_back * np.sin(theta_angle_back * np.pi / 180) * dtheta * 2 * np.pi / 360)
            Q_ext_ADDA[k] = (np.pi * (a / 1000) ** 2 / particle_cross_section) * (4 / x[k] ** 2) * (S11_muel_for[0])  # Compute Qext

        C_sca_for_ADDA[k] = (wavelength[k] / 1000) ** 2 / (2 * np.pi) * np.sum(
            S11_muel_for * np.sin(theta_angle_for * np.pi / 180) * dtheta * 2 * np.pi / 360) / sub_n[k] ** 2  # Integrates forward scattering cross section. Output is in micrometers
        C_sca_back_ADDA[k] = (wavelength[k] / 1000) ** 2 / (2 * np.pi) * np.sum(
            S11_muel_back * np.sin(theta_angle_back * np.pi / 180) * dtheta * 2 * np.pi / 360)  # Integrates backward scattering cross section. Output is in micrometers
        C_ext_ADDA[k] = (4 / x[k] ** 2) * (S11_muel_for[0]) * np.pi * a ** 2 / (1000000)  # Still experimental, it is not working properly...

        C_sca_ADDA_int_tot[k] = C_sca_back_ADDA[k] + C_sca_for_ADDA[k]

        # Mie scattering efficiencies assuming spheres
        if size_type == "Size Along X axis":
            Q_ext_Mie[k] = Qext_Mie(PS_n[k] - 1j * PS_k[k], wavelength[k], a / 2)
            Q_sca_Mie_for[k] = QSca_Mie_for(PS_n[k] - 1j * PS_k[k], wavelength[k], a / 2)
            Q_sca_Mie_back[k] = QSca_Mie_back(PS_n[k] - 1j * PS_k[k], wavelength[k], a / 2)
        else:
            Q_ext_Mie[k] = Qext_Mie(PS_n[k] - 1j * PS_k[k], wavelength[k], a)
            Q_sca_Mie_for[k] = QSca_Mie_for(PS_n[k] - 1j * PS_k[k], wavelength[k], a)
            Q_sca_Mie_back[k] = QSca_Mie_back(PS_n[k] - 1j * PS_k[k], wavelength[k], a)
    print('Upload from ADDA successful')

    # Stacks results and saves them inside SimResfolder

    resultQ = np.column_stack((wavelength, Q_ext_ADDA_auto, Q_abs_ADDA_auto, Q_ext_ADDA, Q_sca_for_ADDA, Q_sca_back_ADDA, Q_sca_ADDA_int, Q_ext_Mie, Q_sca_Mie_for, Q_sca_Mie_back))
    resultC = np.column_stack((wavelength, C_ext_ADDA_auto, C_abs_ADDA_auto, C_sca_ADDA_int, C_sca_ADDA_int_tot, C_sca_for_ADDA, C_sca_back_ADDA))

    filename_Q = f"Simulation_result_Q_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_{Flag}.txt"
    full_pathQ = os.path.join(SimResfolder, filename_Q)
    np.savetxt(full_pathQ, resultQ, fmt='%.6f', header='wavelength Q_ext_ADDA_auto Q_abs_ADDA_auto Q_ext_ADDA Q_sca_for_ADDA Q_sca_back_ADDA Q_sca_ADDA_int Q_ext_Mie Q_sca_Mie_for Q_sca_Mie_back')

    filename_C = f"Simulation_result_C_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_{Flag}.txt"
    full_pathC = os.path.join(SimResfolder, filename_C)
    np.savetxt(full_pathC, resultC, fmt='%.6f', header='wavelength C_ext_ADDA_auto C_abs_ADDA_auto C_sca_ADDA_int C_sca_ADDA_int_tot C_sca_ADDA_int_tot_for C_sca_ADDA_int_tot_back')

    print(f'txt files saved in {SimResfolder}')

    return

def upload_grid(a):  # Upload values from a simulation, calculates scattering quantities and saves them in a .txt file. Called when full phi dependance is required
    # Handle the run parameters
    if wavelength.shape[0] > 1:
        initial_run = int(wavelength[0])
        step_lambda = int(wavelength[1] - wavelength[0])
        final_run = initial_run + step_lambda * wavelength.shape[0]
    else:
        initial_run = int(wavelength[0])
        step_lambda = 1  # no step if only one point
        final_run = initial_run + step_lambda

    S11_muel_grid = np.zeros((len(theta_grid) * len(phi_grid), wavelength.shape[0]))  # Initialize arrays
    S11_theta = np.zeros((len(theta_grid) * len(phi_grid), wavelength.shape[0]))
    S11_phi = np.zeros((len(theta_grid) * len(phi_grid), wavelength.shape[0]))

    if Csca_option == "yes":
        Pol_Y = np.zeros((6, wavelength.shape[0]))
        Pol_X = np.zeros((6, wavelength.shape[0]))
    if Csca_option == "no":
        Pol_Y = np.zeros((4, wavelength.shape[0]))
        Pol_X = np.zeros((4, wavelength.shape[0]))

    Q_ext_ADDA_autoY = np.zeros(wavelength.shape[0])
    Q_abs_ADDA_autoY = np.zeros(wavelength.shape[0])
    Q_ext_ADDA_autoX = np.zeros(wavelength.shape[0])
    Q_abs_ADDA_autoX = np.zeros(wavelength.shape[0])
    Q_sca_ADDA_intY = np.zeros(wavelength.shape[0])
    Q_sca_ADDA_intX = np.zeros(wavelength.shape[0])

    C_ext_ADDA_autoY = np.zeros(wavelength.shape[0])
    C_abs_ADDA_autoY = np.zeros(wavelength.shape[0])
    C_ext_ADDA_autoX = np.zeros(wavelength.shape[0])
    C_abs_ADDA_autoX = np.zeros(wavelength.shape[0])
    C_sca_ADDA_intY = np.zeros(wavelength.shape[0])
    C_sca_ADDA_intX = np.zeros(wavelength.shape[0])

    for j, wl in enumerate(wavelength):
        i = int(wl)  # wavelength value, used for folder naming

        if i < 10:
            run = "00" + str(i)
        if 10 <= i < 100:
            run = "0" + str(i)
        else:
            run = str(i)

            # j = int((i - initial_run)/step_lambda)

        filename_parts = f"run{run}_{particle_shape}_sub{substrate_presence}_g{number_of_dipoles_in_a_row}_m{str(PS_n[j])}_R{str(a)}_{Flag}"
        MuelMatrxGrid = Mainfolder / filename_parts / "mueller_scatgrid"

        S11_muel_grid[:, j] = np.loadtxt(MuelMatrxGrid, delimiter=" ", skiprows=1, usecols=2)
        S11_theta[:, j] = np.loadtxt(MuelMatrxGrid, delimiter=" ", skiprows=1, usecols=0)
        S11_phi[:, j] = np.loadtxt(MuelMatrxGrid, delimiter=" ", skiprows=1, usecols=1)

        CrossSecY = Mainfolder / filename_parts / "CrossSec-Y"  # "_R" + str(a) +

        CrossSecX_candidate = Mainfolder / filename_parts / "CrossSec-X"  # X cross sections

        # If CrossSec-X exists, use it. Otherwise, fall back to CrossSec-Y (for example due to symmetry)
        if os.path.exists(CrossSecX_candidate):
            CrossSecX = CrossSecX_candidate
        else:
            CrossSecX = CrossSecY

        columnY = []
        with open(CrossSecY) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 3:
                    columnY.append([float(parts[2])])

        columnX = []
        with open(CrossSecX) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 3:
                    columnX.append([float(parts[2])])

        Pol_Y[:, j] = np.array(columnY).flatten()
        Pol_X[:, j] = np.array(columnX).flatten()

        Q_abs_ADDA_autoY[j] = Pol_Y[3, j]
        Q_ext_ADDA_autoY[j] = Pol_Y[1, j]
        Q_abs_ADDA_autoX[j] = Pol_X[3, j]
        Q_ext_ADDA_autoX[j] = Pol_X[1, j]

        C_abs_ADDA_autoY[j] = Pol_Y[2, j]
        C_ext_ADDA_autoY[j] = Pol_Y[0, j]
        C_abs_ADDA_autoX[j] = Pol_X[2, j]
        C_ext_ADDA_autoX[j] = Pol_X[0, j]

        if Csca_option == "yes":
            Q_sca_ADDA_intY[j] = Pol_Y[5, j]  # these output are shown only if -Csca option is selected and it is useful to check that everything is correct.
            Q_sca_ADDA_intX[j] = Pol_X[5, j]

            C_sca_ADDA_intY[j] = Pol_Y[4, j]
            C_sca_ADDA_intX[j] = Pol_X[4, j]
        else:
            Q_sca_ADDA_intY[j] = 0
            Q_sca_ADDA_intX[j] = 0

            C_sca_ADDA_intY[j] = 0
            C_sca_ADDA_intX[j] = 0

    Q_abs_ADDA_auto = (Q_abs_ADDA_autoY + Q_abs_ADDA_autoX) / 2
    Q_ext_ADDA_auto = (Q_ext_ADDA_autoY + Q_ext_ADDA_autoX) / 2
    Q_sca_ADDA_int = (Q_sca_ADDA_intY + Q_sca_ADDA_intX) / 2

    C_abs_ADDA_auto = (C_abs_ADDA_autoY + C_abs_ADDA_autoX) / 2
    C_ext_ADDA_auto = (C_ext_ADDA_autoY + C_ext_ADDA_autoX) / 2
    C_sca_ADDA_int = (C_sca_ADDA_intY + C_sca_ADDA_intX) / 2

    # Load Qext directly from CrossSection calculation of output folder, for each wavelength.

    Q_ext_ADDA = np.zeros(len(wavelength))
    Q_sca_for_ADDA = np.zeros(len(wavelength))
    Q_sca_back_ADDA = np.zeros(len(wavelength))
    C_ext_ADDA = np.zeros(len(wavelength))
    C_sca_for_ADDA = np.zeros(len(wavelength))
    C_sca_back_ADDA = np.zeros(len(wavelength))
    Q_sca_Mie_for = np.zeros(len(wavelength))
    Q_sca_Mie_back = np.zeros(len(wavelength))
    Q_ext_Mie = np.zeros(len(wavelength))
    C_sca_ADDA_int_tot = np.zeros(len(wavelength))
    C_sca_ADDA_int_tot_for = np.zeros(len(wavelength))
    C_sca_ADDA_int_tot_back = np.zeros(len(wavelength))
    x = np.zeros(len(wavelength))

    for k in range(len(wavelength)):
        x[k] = 2 * np.pi * a / wavelength[k]  # Compute size parameter

        dtheta = (S11_theta[len(phi_grid), 0] - S11_theta[0, 0]) * np.pi / 180  # Step of theta, as defined in muel_scatgrid
        dphi = (S11_phi[1, 0] - S11_phi[0, 0]) * np.pi / 180  # Step of phi, as defined in muel_scatgrid
        theta_90_index = int(len(phi_grid) * ((len(theta_grid) - 1)) / 2)  # Defines the row number corresponding to theta=90°

        if k == 0:
            print(f'90° is set as row: {theta_90_index}')

        if substrate_presence == "no":  # while working in substrate mode, theta angle must be reversed to assume that positive z axis direction is the impinging light direction...
            C_sca_ADDA_int_tot_for[k] = (((wavelength[k] / 1000) ** 2) / (4 * np.pi ** 2)) * np.sum(
                S11_muel_grid[:theta_90_index, k] * np.sin(S11_theta[:theta_90_index, k] * np.pi / 180) * dtheta * dphi)
            C_sca_ADDA_int_tot_back[k] = (((wavelength[k] / 1000) ** 2) / (4 * np.pi ** 2)) * np.sum(
                S11_muel_grid[theta_90_index:, k] * np.sin(S11_theta[theta_90_index:, k] * np.pi / 180) * dtheta * dphi)
            C_sca_ADDA_int_tot[k] = (((wavelength[k] / 1000) ** 2) / (4 * np.pi ** 2)) * np.sum(
                S11_muel_grid[:, k] * np.sin(S11_theta[:, k] * np.pi / 180) * dtheta * dphi)  # Total cross section, assuming no substrate

        if substrate_presence == "yes":  # while working in substrate mode, theta angle must be reversed to assume that positive z axis direction is the impinging light direction...
            C_sca_ADDA_int_tot_for[k] = (((wavelength[k] / 1000) ** 2) / (4 * np.pi ** 2)) * np.sum(
                S11_muel_grid[theta_90_index:, k] * np.sin(S11_theta[theta_90_index:, k] * np.pi / 180) * dtheta * dphi) / sub_n[k] ** 2
            C_sca_ADDA_int_tot_back[k] = (((wavelength[k] / 1000) ** 2) / (4 * np.pi ** 2)) * np.sum(
                S11_muel_grid[:theta_90_index, k] * np.sin(S11_theta[:theta_90_index, k] * np.pi / 180) * dtheta * dphi)
            C_sca_ADDA_int_tot[k] = C_sca_ADDA_int_tot_back[k] + C_sca_ADDA_int_tot_for[k]

        if size_type == "Size Along X axis":
            Q_ext_Mie[k] = Qext_Mie(PS_n[k] - 1j * PS_k[k], wavelength[k], a / 2)
            Q_sca_Mie_for[k] = QSca_Mie_for(PS_n[k] - 1j * PS_k[k], wavelength[k], a / 2)
            Q_sca_Mie_back[k] = QSca_Mie_back(PS_n[k] - 1j * PS_k[k], wavelength[k], a / 2)
        else:
            Q_ext_Mie[k] = Qext_Mie(PS_n[k] - 1j * PS_k[k], wavelength[k], a)
            Q_sca_Mie_for[k] = QSca_Mie_for(PS_n[k] - 1j * PS_k[k], wavelength[k], a)
            Q_sca_Mie_back[k] = QSca_Mie_back(PS_n[k] - 1j * PS_k[k], wavelength[k], a)

    Q_sca_for_ADDA = C_sca_ADDA_int_tot_for / particle_cross_section
    Q_sca_back_ADDA = C_sca_ADDA_int_tot_back / particle_cross_section

    print('Upload from ADDA successful')

    # Stacks results and saves them inside SimResfolder

    resultQ = np.column_stack((wavelength, Q_ext_ADDA_auto, Q_abs_ADDA_auto, Q_ext_ADDA, Q_sca_for_ADDA, Q_sca_back_ADDA, Q_sca_ADDA_int, Q_ext_Mie, Q_sca_Mie_for, Q_sca_Mie_back))
    resultC = np.column_stack((wavelength, C_ext_ADDA_auto, C_abs_ADDA_auto, C_sca_ADDA_int, C_sca_ADDA_int_tot, C_sca_ADDA_int_tot_for, C_sca_ADDA_int_tot_back))

    filename_Q = f"Simulation_result_Q_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_fromgrid_{Flag}.txt"
    full_pathQ = os.path.join(SimResfolder, filename_Q)
    np.savetxt(full_pathQ, resultQ, fmt='%.6f', header='wavelength Q_ext_ADDA_auto Q_abs_ADDA_auto Q_ext_ADDA Q_sca_for_ADDA Q_sca_back_ADDA Q_sca_ADDA_int Q_ext_Mie Q_sca_Mie_for Q_sca_Mie_back')

    filename_C = f"Simulation_result_C_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_fromgrid_{Flag}.txt"
    full_pathC = os.path.join(SimResfolder, filename_C)
    np.savetxt(full_pathC, resultC, fmt='%.6f', header='wavelength C_ext_ADDA_auto C_abs_ADDA_auto C_sca_ADDA_int C_sca_ADDA_int_tot C_sca_ADDA_int_tot_for C_sca_ADDA_int_tot_back')

    print(f'txt files saved in {SimResfolder}')
    return

def plot(a):
    a_value = str(a)

    data_input_C = f"Simulation_result_C_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_{Flag}.txt"
    full_pathC = os.path.join(SimResfolder, data_input_C)

    data_input_Q = f"Simulation_result_Q_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_{Flag}.txt"
    full_pathQ = os.path.join(SimResfolder, data_input_Q)

    # Load efficiencies
    Q_ext_ADDA_auto = np.loadtxt(full_pathQ, skiprows=1, usecols=1)
    Q_abs_ADDA_auto = np.loadtxt(full_pathQ, skiprows=1, usecols=2)
    Q_sca_for_ADDA = np.loadtxt(full_pathQ, skiprows=1, usecols=4)
    Q_sca_back_ADDA = np.loadtxt(full_pathQ, skiprows=1, usecols=5)
    Q_ext_Mie = np.loadtxt(full_pathQ, skiprows=1, usecols=7)
    Q_sca_Mie_for = np.loadtxt(full_pathQ, skiprows=1, usecols=8)
    Q_sca_Mie_back = np.loadtxt(full_pathQ, skiprows=1, usecols=9)

    Q_ext_ADDA = Q_ext_ADDA_auto  # skip the definition from mueller matrix (#TODO: FIX and use the CrossSec output directly)

    # Load scattering cross sections
    C_ext_ADDA_auto = np.loadtxt(full_pathC, skiprows=1, usecols=1)
    C_abs_ADDA_auto = np.loadtxt(full_pathC, skiprows=1, usecols=2)
    C_sca_ADDA_int = np.loadtxt(full_pathC, skiprows=1, usecols=3)
    C_sca_ADDA_int_tot = np.loadtxt(full_pathC, skiprows=1, usecols=4)
    C_sca_ADDA_int_tot_for = np.loadtxt(full_pathC, skiprows=1, usecols=5)
    C_sca_ADDA_int_tot_back = np.loadtxt(full_pathC, skiprows=1, usecols=6)

    fig, (ax1) = pylab.subplots(1, 1, figsize=(4, 8), dpi=600)

    ax1.plot(wavelength, Q_ext_Mie, '--', color=(0, 1, 0), label='ext Mie')
    ax1.plot(wavelength, Q_sca_Mie_for, '--', color=(1, 0, 0), label='sca for Mie')
    ax1.plot(wavelength, Q_sca_Mie_back, '--', color=(0, 0.3, 1), label='sca back Mie')
    ax1.plot(wavelength, Q_ext_ADDA, color='green', label='ext ADDA auto')
    ax1.plot(wavelength, Q_sca_for_ADDA, color=(0.5, 0, 0), label='sca for ADDA')
    ax1.plot(wavelength, Q_sca_back_ADDA, color=(0, 0, 0.5), label='sca back ADDA')
    ax1.legend(loc="upper right", fontsize="small", frameon=True, framealpha=0)
    ax1.set_title(f"Efficiency factors ADDA vs Mie, R={a_value}nm")

    fig2, (ax4) = pylab.subplots(1, 1, figsize=(4, 8), dpi=600)
    ax4.plot(wavelength, C_sca_ADDA_int_tot_for, color=(0, 0, 0), label='Csca for')
    ax4.plot(wavelength, C_sca_ADDA_int_tot_back, color=(1, 0, 0), label='Csca back')
    ax4.plot(wavelength, C_sca_ADDA_int_tot, color=(0, 0.8, 0), label='Csca for+Csca back ')
    ax4.plot(wavelength, C_sca_ADDA_int, '--', color=(0, 0, 1), label='Csca out ADDA')
    ax4.plot(wavelength, C_abs_ADDA_auto, color=(1, 0, 1), label='Cabs auto')
    ax4.plot(wavelength, C_ext_ADDA_auto, color=(0, 0.8, 1), label='Cext auto')
    ax4.legend(loc="upper right", fontsize="small", frameon=True, framealpha=0)
    ax4.set_title(f"ADDA Scattering cross sections, R={a_value}nm")

    pylab.show()

def plot_grid(a):
    a_value = str(a)

    data_input_C = f"Simulation_result_C_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_fromgrid_{Flag}.txt"
    full_pathC = os.path.join(SimResfolder, data_input_C)

    data_input_Q = f"Simulation_result_Q_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_fromgrid_{Flag}.txt"
    full_pathQ = os.path.join(SimResfolder, data_input_Q)

    # Load efficiencies
    Q_ext_ADDA_auto = np.loadtxt(full_pathQ, skiprows=1, usecols=1)
    Q_abs_ADDA_auto = np.loadtxt(full_pathQ, skiprows=1, usecols=2)
    Q_sca_for_ADDA = np.loadtxt(full_pathQ, skiprows=1, usecols=4)
    Q_sca_back_ADDA = np.loadtxt(full_pathQ, skiprows=1, usecols=5)
    Q_ext_Mie = np.loadtxt(full_pathQ, skiprows=1, usecols=7)
    Q_sca_Mie_for = np.loadtxt(full_pathQ, skiprows=1, usecols=8)
    Q_sca_Mie_back = np.loadtxt(full_pathQ, skiprows=1, usecols=9)

    Q_ext_ADDA = Q_ext_ADDA_auto  # skip the definition from mueller matrix (#TODO: FIX) and use the CrossSec output directly

    # Load scattering cross sections
    C_ext_ADDA_auto = np.loadtxt(full_pathC, skiprows=1, usecols=1)
    C_abs_ADDA_auto = np.loadtxt(full_pathC, skiprows=1, usecols=2)
    C_sca_ADDA_int = np.loadtxt(full_pathC, skiprows=1, usecols=3)
    C_sca_ADDA_int_tot = np.loadtxt(full_pathC, skiprows=1, usecols=4)
    C_sca_ADDA_int_tot_for = np.loadtxt(full_pathC, skiprows=1, usecols=5)
    C_sca_ADDA_int_tot_back = np.loadtxt(full_pathC, skiprows=1, usecols=6)

    fig, (ax1) = pylab.subplots(1, 1, figsize=(4, 6), dpi=600)

    ax1.plot(wavelength, Q_ext_Mie, '--', color=(0, 1, 0), label='ext Mie')
    ax1.plot(wavelength, Q_sca_Mie_for, '--', color=(1, 0, 0), label='sca for Mie')
    ax1.plot(wavelength, Q_sca_Mie_back, '--', color=(0, 0.3, 1), label='sca back Mie')
    ax1.plot(wavelength, Q_ext_ADDA, color='green', label='ext ADDA auto')
    ax1.plot(wavelength, Q_sca_for_ADDA, color=(0.5, 0, 0), label='sca for ADDA')
    ax1.plot(wavelength, Q_sca_back_ADDA, color=(0, 0, 0.5), label='sca back ADDA')
    ax1.set_xlabel("Wavelength (nm)", fontsize=12, fontweight="bold")
    ax1.set_ylabel("Efficiency factor", fontsize=12, fontweight="bold")
    ax1.legend(loc="upper right", fontsize="small", frameon=True, framealpha=0)
    ax1.set_title(f"Efficiency factors ADDA vs Mie, R={a_value}nm")

    fig2, (ax4) = pylab.subplots(1, 1, figsize=(4, 6), dpi=600)
    ax4.plot(wavelength, C_sca_ADDA_int_tot_for, color=(0, 0, 0), label='Csca for')
    ax4.plot(wavelength, C_sca_ADDA_int_tot_back, color=(1, 0, 0), label='Csca back')
    ax4.plot(wavelength, C_sca_ADDA_int_tot, color=(0, 0.8, 0), label='Csca for+Csca back ')
    ax4.plot(wavelength, C_sca_ADDA_int, '--', color=(0, 0, 1), label='Csca out ADDA')
    ax4.plot(wavelength, C_abs_ADDA_auto, color=(1, 0, 1), label='Cabs auto')
    ax4.plot(wavelength, C_ext_ADDA_auto, color=(0, 0.8, 1), label='Cext auto')
    ax4.legend(loc="upper right", fontsize="small", frameon=True, framealpha=0)
    ax4.set_xlabel("Wavelength (nm)", fontsize=12, fontweight="bold")
    ax4.set_ylabel("Cross section (µm²)", fontsize=12, fontweight="bold")
    ax4.set_title(f"ADDA Scattering cross sections, R={a_value}nm")

    pylab.show()

# --- 6. GUI FUNCTIONS ---

def load_input_data():
    global wavelength, PS_n, PS_k, sub_n
    filename = entry_input_file.get()
    try:
        if not os.path.isabs(filename): filename = Mainfolder / filename
        data = np.loadtxt(filename)
        data = np.atleast_2d(data)
        wavelength = data[:, 0]; PS_n = data[:, 1]; PS_k = data[:, 2]; sub_n = data[:, 3]
        return True
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load input file: {e}")
        return False

def browse_input_file():
    f = filedialog.askopenfilename(initialdir=Mainfolder)
    if f: entry_input_file.delete(0, tk.END); entry_input_file.insert(0, f)

def get_gui_params():
    global phi_stop, theta_stop, keep_bat_files, size_type, phi_grid, particle_cross_section
    global store_scat_grid_option, Csca_option, particle_shape, substrate_presence
    global grid_value, number_of_dipoles_in_a_row, theta_grid, parallel_processes
    global ntheta_value, custom_shape, distance_from_surf, Flag
    global GPU_acc
    
    grid_value = int(entry_grid.get())
    number_of_dipoles_in_a_row = str(grid_value)
    substrate_presence = substrate_var.get()
    a = float(entry_a.get())
    particle_shape = entry_shape.get()
    parallel_processes = int(entry_parallel.get())
    
    theta_step = float(entry_theta_step.get())
    theta_stop = float(entry_theta_stop.get())
    theta_grid = np.arange(0, theta_stop + theta_step/2, theta_step)
    ntheta_value = len(theta_grid) - 1
    
    phi_step = float(entry_phi_step.get())
    phi_stop = float(entry_phi_stop.get())
    phi_grid = np.arange(0, phi_stop + phi_step/2, phi_step)
    
    particle_cross_section = np.pi * (a/1000)**2
    
    GPU_acc=GPU_acc_var.get()
    custom_shape = custom_shape_var.get()
    distance_from_surf = distance_from_surf_var.get()
    Flag = entry_Flag.get()
    Csca_option = Csca_option_var.get()
    store_scat_grid_option = store_scat_grid_option_var.get()
    size_type = size_type_var.get()
    keep_bat_files = keep_bat_files_var.get()
    
    return a

def run_upload_plot_sim():
    if not load_input_data(): return
    if not setup_adda_environment(): return
    try:
        a = get_gui_params()
        if store_scat_grid_option == "yes":
            go(a); upload_grid(a); plot_grid(a)
        else:
            go(a); upload(a); plot(a)
        messagebox.showinfo("Done", "Simulation finished.")
    except Exception as e: messagebox.showerror("Error", str(e))

def run_sim():
    if not load_input_data() or not setup_adda_environment(): return
    try:
        a = get_gui_params()
        go(a)
        messagebox.showinfo("Done", "Run finished.")
    except Exception as e: messagebox.showerror("Error", str(e))

def upload_sim():
    if not load_input_data(): return
    try:
        a = get_gui_params()
        if store_scat_grid_option_var.get() == "yes": upload_grid(a)
        else: upload(a)
        messagebox.showinfo("Done", "Upload finished.")
    except Exception as e: messagebox.showerror("Error", str(e))

def plot_sim():
    if not load_input_data(): return
    try:
        a = get_gui_params()
        if store_scat_grid_option_var.get() == "yes": plot_grid(a)
        else: plot(a)
    except Exception as e: messagebox.showerror("Error", str(e))

def kill_folder_gui():
    if not load_input_data(): return
    try:
        a = get_gui_params()
        kill_folder(a)
        messagebox.showinfo("Done", "Folders deleted.")
    except Exception as e: messagebox.showerror("Error", str(e))


# --- 7. MAIN GUI LAYOUT ---

root = tk.Tk()
root.title("ADDA Simulation GUI v5.0")
row = 0

# 1. Input File
tk.Label(root, text="Input File:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
entry_input_file = tk.Entry(root, width=40); entry_input_file.insert(0, "lambda_nk_test.txt")
entry_input_file.grid(row=row, column=1, columnspan=2, sticky="ew", pady=2)
tk.Button(root, text="Browse...", command=browse_input_file).grid(row=row, column=3, sticky="w", padx=5); row+=1

# 2. Parallel
tk.Label(root, text="Parallel Processes:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
entry_parallel = tk.Entry(root, width=6); entry_parallel.insert(0, "4")
entry_parallel.grid(row=row, column=1, sticky="w"); row+=1

# 3. Grid Value
tk.Label(root, text="Grid Value:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
entry_grid = tk.Entry(root, width=6); entry_grid.insert(0, "32")
entry_grid.grid(row=row, column=1, sticky="w"); row+=1

# 4. Particle Size
tk.Label(root, text="Particle dimension [nm]:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
entry_a = tk.Entry(root, width=6); entry_a.insert(0, "50")
entry_a.grid(row=row, column=1, sticky="w"); row+=1

# 5. Dimension Type
tk.Label(root, text="Particle dimension type:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
size_type_var = tk.StringVar(value="Equivalent radius")
ttk.Combobox(root, textvariable=size_type_var, values=["Equivalent radius", "Size Along X axis"], width=20).grid(row=row, column=1, columnspan=2, sticky="w"); row+=1

# 6. Substrate
tk.Label(root, text="Substrate:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
substrate_var = tk.StringVar(value="no")
ttk.Combobox(root, textvariable=substrate_var, values=["yes", "no"], width=5).grid(row=row, column=1, sticky="w"); row+=1

# 7. Distance
tk.Label(root, text="Distance from surface [nm]:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
distance_from_surf_var = tk.Entry(root, width=6); distance_from_surf_var.insert(0, "50")
distance_from_surf_var.grid(row=row, column=1, sticky="w"); row+=1

# 8. Custom Shape
tk.Label(root, text="Custom Shape:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
custom_shape_var = tk.StringVar(value="no")
ttk.Combobox(root, textvariable=custom_shape_var, values=["yes", "no"], width=5).grid(row=row, column=1, sticky="w"); row+=1

# 9. Particle Shape (Standard)
tk.Label(root, text="Particle Shape:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
entry_shape = tk.Entry(root, width=20); entry_shape.insert(0, "sphere")
entry_shape.grid(row=row, column=1, columnspan=2, sticky="w"); row+=1

# 10. Theta Stop
tk.Label(root, text="Theta Stop [deg]:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
entry_theta_stop = tk.Entry(root, width=6); entry_theta_stop.insert(0, "180")
entry_theta_stop.grid(row=row, column=1, sticky="w"); row+=1

# 11. Theta Step
tk.Label(root, text="Theta Step [deg]:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
entry_theta_step = tk.Entry(root, width=6); entry_theta_step.insert(0, "0.5")
entry_theta_step.grid(row=row, column=1, sticky="w"); row+=1

# 12. Phi Stop
tk.Label(root, text="Phi Stop [deg]:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
entry_phi_stop = tk.Entry(root, width=6); entry_phi_stop.insert(0, "360")
entry_phi_stop.grid(row=row, column=1, sticky="w"); row+=1

# 13. Phi Step
tk.Label(root, text="Phi Step [deg]:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
entry_phi_step = tk.Entry(root, width=6); entry_phi_step.insert(0, "2")
entry_phi_step.grid(row=row, column=1, sticky="w"); row+=1

# 14. Flag
tk.Label(root, text="Flag:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
entry_Flag = tk.Entry(root, width=20); entry_Flag.insert(0, "")
entry_Flag.grid(row=row, column=1, columnspan=2, sticky="w"); row+=1

# 15. Csca Option
tk.Label(root, text="Internal Csca integration:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
Csca_option_var = tk.StringVar(value="yes")
ttk.Combobox(root, textvariable=Csca_option_var, values=["yes", "no"], width=5).grid(row=row, column=1, sticky="w"); row+=1

# 16. Store Scat Grid
tk.Label(root, text="Enable Phi scattering grid:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
store_scat_grid_option_var = tk.StringVar(value="yes")
ttk.Combobox(root, textvariable=store_scat_grid_option_var, values=["yes", "no"], width=5).grid(row=row, column=1, sticky="w"); row+=1

# 17. Keep Bat
tk.Label(root, text="Keep .bat files:").grid(row=row, column=0, sticky="w", padx=5, pady=2)
keep_bat_files_var = tk.StringVar(value="no")
ttk.Combobox(root, textvariable=keep_bat_files_var, values=["yes", "no"], width=5).grid(row=row, column=1, sticky="w"); row+=1

# 18. GPU OCL acceleration
tk.Label(root, text="GPU accelleration: no, or GPU number").grid(row=row, column=0, sticky="w", padx=5, pady=2)
GPU_acc_var = tk.Entry(root, width=6); GPU_acc_var.insert(0, "no")
GPU_acc_var.grid(row=row, column=1, sticky="w"); row+=1

# Buttons
btn_f = tk.Frame(root)
btn_f.grid(row=row, column=0, columnspan=4, pady=10)
tk.Button(btn_f, text="Run, Upload & Plot", command=run_upload_plot_sim, width=20, bg="#ddffdd").pack(pady=2)
tk.Button(btn_f, text="Run", command=run_sim, width=20).pack(pady=2)
tk.Button(btn_f, text="Upload", command=upload_sim, width=20).pack(pady=2)
tk.Button(btn_f, text="Plot", command=plot_sim, width=20).pack(pady=2)
tk.Button(btn_f, text="Kill Folders", command=kill_folder_gui, width=20, bg="#ffdddd").pack(pady=5)

# --- 8. STARTUP ---
setup_adda_environment() # Auto-detect adda.exe
root.mainloop()