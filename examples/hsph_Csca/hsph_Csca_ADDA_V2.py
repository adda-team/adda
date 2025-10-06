"""
V2.0

@author: Mattia Andrini, mattia.andrini@unicatt.it

You can write me for anything about the code. Suggestions are welcome.

I do not claim myself to be a good programmer. Please be kind :D

A paper making use of this version of the code is currently in preparation. 
In the meantime, if you use this code, please cite (https://doi.org/10.1021/acs.analchem.5c00168), as its core functions were originally developed there.

"""
import time
import numpy as np
import pylab
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import tkinter as tk
from tkinter import ttk, messagebox
import shutil

'This code allows to run controlled simualtion with the software ADDA given a wavelength spectrum input, and integrate on hemispherical angle domains.'

'First, set the Mainfolder and SimResfolder, and put this script, the input_filename and all_dir_params.dat in the Mainfolder.' 
'Then, run this script. A GUI window should appear.'
'Clik on Run to start a simulation, on Upload to perform the integration and save the results, on Plot to visualize the results.'
'Read the guide for more details.'

Mainfolder="your path to\\adda-1.4.0\\win64\\" #set name of the working folder
SimResfolder = r"your path to\\SimulationResults" #specify where to save the scattering results

### PATH EXAMPLE ###
#Mainfolder="C:\\Users\\m51vo\\Desktop\\adda-1.4.0\\win64\\"
#SimResfolder = "C:\\Users\\m51vo\\Desktop\\adda-1.4.0\\SimulationResults"
### END PATH EXAMPLE ###

input_filename="lambda_nk_test.txt"   #4 columns of values separated by spaces in this order: wavelength, real part of the scatterer refractive index (n), imaginary part of the scatterer refractive index (k), real part of the substrate refractive index (even if you are working in free space mode)"
                          
data = np.loadtxt(input_filename)                   #Load the data from the file
data = np.atleast_2d(data)                          #Make sure it's at least 2D so slicing works even for a single row

wavelength = data[:, 0]   # wavelength values
PS_n       = data[:, 1]   # real part of RI
PS_k       = data[:, 2]   # imaginary part of RI
sub_n      = data[:, 3]   # substrate RI


def go(a):                  #main function to perform a simulation for a given size "a" (radius if spheres, otherwise check the actual dimension) for a full spectrum specified in the .txt file "input_filename"
    start_time=time.time()                  #needed to print the time needed for the simulation
    bat_files = []
    for k in range(wavelength.shape[0]):    #for each wavelength in the input_file, it creates a batch file with the lambda and refractive index of the k-row corresponding to that lambda.
        bat_file=generate_bat_file(wavelength[k],PS_n[k],PS_k[k],sub_n[k],k,a)
        bat_files.append((bat_file, k))
                                            #batch files are run one each time, for each wavelength.
     # 2. Run batch files in parallel
    print("Starting parallel simulations...") #used to perform simulations in parallel, the time used for each wavelength is printed in the console.
    with ThreadPoolExecutor(max_workers=parallel_processes) as executor:
        futures = [executor.submit(run_batch_file, bat_file, k) for bat_file, k in bat_files]
        for future in as_completed(futures):
            future.result()    
    if keep_bat_files=="no":
        for bat_file, _ in bat_files:  #removes bat files at the end of the simulation. comment this "for" if you have troubles to keep the .bat files for debugging. add "pause" at the end of the .bat file to read the .cmd comments and see what went wrong.
            try:
                os.remove(bat_file)
            except Exception as e:  
                print(f"Could not delete {bat_file}: {e}")
    print(f"R={a}nm simulation took: {time.time() - start_time:.2f} seconds")
    print('Simulation performed!')
    return
def generate_bat_file(wavelength,PS_n,PS_k,sub_n,k,a):     #first variable is the wavelength which is converted to micrometers (ADDA takes µm). The output folder will display wavelenth, particle shape, substrate presence, grid, particle n, dimension and an user-dependant flag to assure unicity.
    bat_filename = f"run_{k}.bat"    
    with open(bat_filename, 'w') as f:
        # --- Header --- set working folder and avoid printing the log.
        f.write("echo off\n")
        f.write("cd win64\n")
        
       # --- Common parameters ---
        wavelength_micro = float(wavelength)/1000   # nm → µm
        lambda_value = int(wavelength)
        eq_rad_value = float(a)/1000
        distance_from_surf_value = float(distance_from_surf)/1000
        command =  f"adda -scat_matr muel -lambda {wavelength_micro} -m {PS_n} {PS_k} "  #base command. Use adda_ocl -gpu [value] to exploit double precision if you have a dedicated GPU.
        
        if size_type=="Equivalent radius":
            command += f"-eq_rad {eq_rad_value} "
        if size_type=="Size Along X axis":
            command += f"-size {eq_rad_value} "    
        
        if substrate_presence == "yes":    #surface mode
            command += f"-surf {distance_from_surf_value} {sub_n} 0 -prop 0 0 -1 " #surface mode.
        if substrate_presence == "no":
            command += "-prop 0 0 1 "
       
        if custom_shape == "no":
            command += f"-shape {particle_shape} -grid {grid_value} " #you may need to specify more parameters if you are not working with spheres. For example if particle_shape=ellipsoid" --> -shape {particle_shape} [y/x] [z/x] "
        if custom_shape == "yes":
            command += f"-shape read {particle_shape}.geom "   
       
        if Csca_option== "yes": #enable ADDA internal integration
            command+="-Csca "
              
        if store_scat_grid_option=="yes":
            command += "-store_scat_grid "
       
        if store_scat_grid_option=="no":
            command += f"-ntheta {ntheta_value} "
               
        command += f"-dir run{lambda_value}_{particle_shape}_sub{substrate_presence}_g{grid_value}_m{PS_n}_R{str(a)}_{Flag} \n " #set the name of the folder. if you change this, you will need to update the folder name in upload function.                    
        f.write(command)
            
    return bat_filename    
    
    
def run_batch_file(bat_file, k):      #just runs the .bat file and prints the simulation time for each wavelength. not intended to be used alone. 
    start_time = time.time()
    try:            #to print time needed to simulate each batch file. redundant?
        result = subprocess.run(bat_file, shell=True, check=True)
        end_time = time.time()
        elapsed = end_time - start_time
        #print(f"Batch file completed with return code: {result.returncode}")
        print(f"{wavelength[k]} simulation took: {elapsed:.2f} seconds")
    except subprocess.CalledProcessError as e:
        print(f"Error running batch file: {e}")


def upload(a):  #upload values from a simulation, calculates scattering quantities and saves them in a .txt file when the phi scatering option is disabled. 
    
# Handle the run parameters
    if wavelength.shape[0] > 1:
        initial_run = int(wavelength[0])
        step_lambda = int(wavelength[1] - wavelength[0])
        final_run = initial_run + step_lambda * wavelength.shape[0]
    else:
        initial_run = int(wavelength[0])
        step_lambda = 1   # no step if only one point
        final_run = initial_run+step_lambda     #function discontinued
    
    S11_muel = np.zeros((len(theta_grid), wavelength.shape[0]))   #just inizialize the arrays
    
    if Csca_option=="yes":
        Pol_Y = np.zeros((6, wavelength.shape[0]))
        Pol_X = np.zeros((6, wavelength.shape[0]))
    if Csca_option=="no":
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
        
        if i<10:
            run = "00" + str(i)
        if i>=10 and i<100:
            run = "0" + str(i)
        else:
            run = str(i)       
        
        MuelMatrx = Mainfolder + "run" + run + "_" + \
        particle_shape + "_sub"+ substrate_presence+ "_g" + number_of_dipoles_in_a_row + \
        "_m" + str(PS_n[j]) + "_R" + str(a)+ f"_{Flag}" +"\\mueller"   #location of mueller matrix
        
        
        #j = int((i - initial_run)/step_lambda)
        
        
        if wavelength.shape[0]==1:
            theta = np.loadtxt(MuelMatrx, delimiter=" ", skiprows=1, usecols=0)[:len(theta_grid)]
                
        if i == initial_run:
            theta = np.loadtxt(MuelMatrx, delimiter=" ", skiprows=1, usecols=0)[:len(theta_grid)]
        S11_muel[:len(theta),j] = np.loadtxt(MuelMatrx, delimiter=" ", skiprows=1, usecols=1)[:len(theta)]
        
        CrossSecY = Mainfolder + "run" +run + "_" + \
        particle_shape + "_sub"+ substrate_presence+ "_g" + number_of_dipoles_in_a_row + \
        "_m" + str(PS_n[j])+ "_R" + str(a)+ f"_{Flag}" +"\\CrossSec-Y"         #Y cross sections
              
        CrossSecX_candidate = Mainfolder + "run" + run + "_" + \
        particle_shape + "_sub"+ substrate_presence+ "_g" + number_of_dipoles_in_a_row + \
        "_m" + str(PS_n[j])+ "_R" + str(a)+ f"_{Flag}" +"\\CrossSec-X"    #X cross sections

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
        
        Pol_Y[:,j] = np.array(columnY).flatten()
        Pol_X[:,j] = np.array(columnX).flatten()
        
        #loads cross sections and efficiency factors from the CrossSec files, for each wavelength index j. 
        
        Q_abs_ADDA_autoY[j]=Pol_Y[3,j]
        Q_ext_ADDA_autoY[j]=Pol_Y[1,j]
        Q_abs_ADDA_autoX[j]=Pol_X[3,j] 
        Q_ext_ADDA_autoX[j]=Pol_X[1,j]
        
        C_abs_ADDA_autoY[j]=Pol_Y[2,j] 
        C_ext_ADDA_autoY[j]=Pol_Y[0,j]
        C_abs_ADDA_autoX[j]=Pol_X[2,j] 
        C_ext_ADDA_autoX[j]=Pol_X[0,j]
        
        if Csca_option=="yes":        
            Q_sca_ADDA_intY[j]=Pol_Y[5,j] #these output are shown only if -Csca option is selected and it useful to check that everything is correct.
            Q_sca_ADDA_intX[j]=Pol_X[5,j]
        
            C_sca_ADDA_intY[j]=Pol_Y[4,j] 
            C_sca_ADDA_intX[j]=Pol_X[4,j]
        else:
            Q_sca_ADDA_intY[j]=0
            Q_sca_ADDA_intX[j]=0
        
            C_sca_ADDA_intY[j]=0 
            C_sca_ADDA_intX[j]=0
                       
    Q_abs_ADDA_auto=(Q_abs_ADDA_autoY+Q_abs_ADDA_autoX)/2    #averages for different directions
    Q_ext_ADDA_auto=(Q_ext_ADDA_autoY+Q_ext_ADDA_autoX)/2    
    Q_sca_ADDA_int=(Q_sca_ADDA_intY+Q_sca_ADDA_intX)/2
    
    C_abs_ADDA_auto=(C_abs_ADDA_autoY+C_abs_ADDA_autoX)/2
    C_ext_ADDA_auto=(C_ext_ADDA_autoY+C_ext_ADDA_autoX)/2    
    C_sca_ADDA_int=(C_sca_ADDA_intY+C_sca_ADDA_intX)/2
        
    theta_90_lab=int(np.round(len(theta_grid)/2)+1)
    theta_91_lab_lab=int(np.round(len(theta_grid)/2)+2)
    
    print(f'90° is set as row {theta_90_lab}')         #prints the row of the mueller file assigned to theta=90°. useful to check consistency.
    
    if substrate_presence=="no":
        theta_angle_for = theta[:theta_90_lab] #defines set of angle for forward scattering (0-90°). 361 label correspond to theta=90° if ngrid=720 (step=0.25°). otherwise change
        theta_angle_back = theta[theta_91_lab_lab:]
    if substrate_presence=="yes": #while working in substrate mode, theta angle must be reversed to assume that positive z axis direction is the impinging light direction...
        theta_angle_for = theta[theta_91_lab_lab:] 
        theta_angle_back = theta[:theta_90_lab]
            
    Q_ext_ADDA= np.zeros(len(wavelength))
    Q_sca_for_ADDA=np.zeros(len(wavelength))
    Q_sca_back_ADDA=np.zeros(len(wavelength))
    C_ext_ADDA= np.zeros(len(wavelength))    
    C_sca_for_ADDA=np.zeros(len(wavelength))
    C_sca_back_ADDA=np.zeros(len(wavelength))
    Q_sca_Mie_for=np.zeros(len(wavelength))
    Q_sca_Mie_back=np.zeros(len(wavelength))
    Q_ext_Mie=np.zeros(len(wavelength))
    C_sca_ADDA_int_tot=np.zeros(len(wavelength))

    x=np.zeros(len(wavelength))
    
    for k in range(len(wavelength)):
        x[k]=2*np.pi*a/wavelength[k]    #compute size parameter
        
        if substrate_presence=="no":
            S11_muel_for = S11_muel[:theta_90_lab, k]
            S11_muel_back = S11_muel[theta_91_lab_lab:, k]
            
        if substrate_presence=="yes": #while working in substrate mode, theta angle must be reversed to assume that positive z axis direction is the impinging light direction...
            S11_muel_for = S11_muel[theta_91_lab_lab:, k]
            S11_muel_back = S11_muel[:theta_90_lab:, k]            
        
        if substrate_presence=="no":  #in substrate mode, the scattering effeciency of the light propagating in the substrate must be divided by the refr index of the subatrate squared.
            sub_n[k]=1
        
        dtheta=theta_grid[1]-theta_grid[0] #takes the theta step from the mueller file.
        
        if size_type=="Size Along X axis":
            Q_sca_for_ADDA[k] = (np.pi*(a/1000)**2/particle_cross_section)*4*2*x[k]**(-2)*np.sum(S11_muel_for*np.sin(theta_angle_for*np.pi/180)*dtheta*2*np.pi/360)/sub_n[k]**2    #integrates forward scattering efficencies assuming SPHERES!
            Q_sca_back_ADDA[k] = (np.pi*(a/1000)**2/particle_cross_section)*4*2*x[k]**(-2)*np.sum(S11_muel_back*np.sin(theta_angle_back*np.pi/180)*dtheta*2*np.pi/360)
            Q_ext_ADDA[k]=(np.pi*(a/1000)**2/particle_cross_section)*4*(4/x[k]**2)*(S11_muel_for[0]) #compute Qext
        else:
            Q_sca_for_ADDA[k] = (np.pi*(a/1000)**2/particle_cross_section)*2*x[k]**(-2)*np.sum(S11_muel_for*np.sin(theta_angle_for*np.pi/180)*dtheta*2*np.pi/360)/sub_n[k]**2    #integrates forward scattering efficencies assuming SPHERES!
            Q_sca_back_ADDA[k] = (np.pi*(a/1000)**2/particle_cross_section)*2*x[k]**(-2)*np.sum(S11_muel_back*np.sin(theta_angle_back*np.pi/180)*dtheta*2*np.pi/360)
            Q_ext_ADDA[k]=(np.pi*(a/1000)**2/particle_cross_section)*(4/x[k]**2)*(S11_muel_for[0]) #compute Qext
        
        
        C_sca_for_ADDA[k] = (wavelength[k]/1000)**2/(2*np.pi)*np.sum(S11_muel_for*np.sin(theta_angle_for*np.pi/180)*dtheta*2*np.pi/360)/sub_n[k]**2 #integrates forward scattering cross section. output is in micrometers
        C_sca_back_ADDA[k] = (wavelength[k]/1000)**2/(2*np.pi)*np.sum(S11_muel_back*np.sin(theta_angle_back*np.pi/180)*dtheta*2*np.pi/360) #integrates backward scattering cross section. output is in micrometers
        C_ext_ADDA[k]=(4/x[k]**2)*(S11_muel_for[0])*np.pi*a**2/(1000000) #still experimental, is not working properly...
        
        C_sca_ADDA_int_tot[k]=C_sca_back_ADDA[k]+C_sca_for_ADDA[k]
        
        #Mie scattering efficiencies assuming spheres
        if size_type=="Size Along X axis":
            Q_ext_Mie[k]=Qext_Mie(PS_n[k]-1j*PS_k[k],wavelength[k],a/2)
            Q_sca_Mie_for[k]=QSca_Mie_for(PS_n[k]-1j*PS_k[k],wavelength[k],a/2)
            Q_sca_Mie_back[k]=QSca_Mie_back(PS_n[k]-1j*PS_k[k],wavelength[k],a/2)
        else:
            Q_ext_Mie[k]=Qext_Mie(PS_n[k]-1j*PS_k[k],wavelength[k],a)
            Q_sca_Mie_for[k]=QSca_Mie_for(PS_n[k]-1j*PS_k[k],wavelength[k],a)
            Q_sca_Mie_back[k]=QSca_Mie_back(PS_n[k]-1j*PS_k[k],wavelength[k],a)            
    print('Upload from ADDA successfull')
    
    #stacks results and save them inside SimResfolder
    
    resultQ=np.column_stack((wavelength,Q_ext_ADDA_auto,Q_abs_ADDA_auto,Q_ext_ADDA,Q_sca_for_ADDA,Q_sca_back_ADDA,Q_sca_ADDA_int,Q_ext_Mie,Q_sca_Mie_for,Q_sca_Mie_back))
    resultC=np.column_stack((wavelength,C_ext_ADDA_auto,C_abs_ADDA_auto,C_sca_ADDA_int,C_sca_ADDA_int_tot,C_sca_for_ADDA,C_sca_back_ADDA))    
    
    filename_Q = f"Simulation_result_Q_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_{Flag}.txt"
    full_pathQ = os.path.join(SimResfolder, filename_Q)
    np.savetxt(full_pathQ,resultQ,fmt='%.6f',header='wavelength Q_ext_ADDA_auto Q_abs_ADDA_auto Q_ext_ADDA Q_sca_for_ADDA Q_sca_back_ADDA Q_sca_ADDA_int Q_ext_Mie Q_sca_Mie_for Q_sca_Mie_back')
    
    filename_C = f"Simulation_result_C_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_{Flag}.txt"
    full_pathC = os.path.join(SimResfolder, filename_C)
    np.savetxt(full_pathC,resultC,fmt='%.6f',header='wavelength C_ext_ADDA_auto C_abs_ADDA_auto C_sca_ADDA_int C_sca_ADDA_int_tot C_sca_ADDA_int_tot_for C_sca_ADDA_int_tot_back')
  
    print(f'txt files saved in {SimResfolder}')  
  
    return

def upload_grid(a):      #upload values from a simulation, calculates scattering quantities and saves them in a .txt file. USE WITH -store_scat_grid
# Handle the run parameters
    if wavelength.shape[0] > 1:
        initial_run = int(wavelength[0])
        step_lambda = int(wavelength[1] - wavelength[0])
        final_run = initial_run + step_lambda * wavelength.shape[0]
    else:
        initial_run = int(wavelength[0])
        step_lambda = 1   # no step if only one point
        final_run = initial_run+step_lambda

    S11_muel_grid = np.zeros((len(theta_grid)*len(phi_grid), wavelength.shape[0]))  #inizialize arrays
    S11_theta = np.zeros((len(theta_grid)*len(phi_grid), wavelength.shape[0]))
    S11_phi = np.zeros((len(theta_grid)*len(phi_grid), wavelength.shape[0]))
    
    if Csca_option=="yes":
        Pol_Y = np.zeros((6, wavelength.shape[0]))
        Pol_X = np.zeros((6, wavelength.shape[0]))
    if Csca_option=="no":
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
        
        if i<10:
            run = "00" + str(i)
        if i>=10 and i<100:
            run = "0" + str(i)
        else:
            run = str(i)       
               
        #j = int((i - initial_run)/step_lambda)
        
        MuelMatrxGrid = Mainfolder +"run" + run + "_" + \
        particle_shape + "_sub"+ substrate_presence+ "_g" + number_of_dipoles_in_a_row + \
        "_m" + str(PS_n[j]) + "_R" + str(a)+ f"_{Flag}" +"\\mueller_scatgrid"

        
        S11_muel_grid[:,j]=np.loadtxt(MuelMatrxGrid,delimiter=" ", skiprows=1, usecols=2)     
        S11_theta[:,j]=np.loadtxt(MuelMatrxGrid,delimiter=" ", skiprows=1, usecols=0)
        S11_phi[:,j]=np.loadtxt(MuelMatrxGrid,delimiter=" ", skiprows=1, usecols=1)
        
        CrossSecY = Mainfolder +"run" + run + "_" + \
        particle_shape + "_sub"+ substrate_presence+ "_g" + number_of_dipoles_in_a_row + \
        "_m" + str(PS_n[j])+ "_R" + str(a)+ f"_{Flag}" +"\\CrossSec-Y" #"_R" + str(a) + 
        
        CrossSecX_candidate = Mainfolder +"run" + run + "_" + \
        particle_shape + "_sub"+ substrate_presence+ "_g" + number_of_dipoles_in_a_row + \
        "_m" + str(PS_n[j])+ "_R" + str(a)+ f"_{Flag}" +"\\CrossSec-X"    #X cross sections

        # If CrossSec-X exists, use it. Otherwise, fall back to CrossSec-Y (example due to symmetry)
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
        
        Pol_Y[:,j] = np.array(columnY).flatten()
        Pol_X[:,j] = np.array(columnX).flatten()
        
        Q_abs_ADDA_autoY[j]=Pol_Y[3,j] 
        Q_ext_ADDA_autoY[j]=Pol_Y[1,j]
        Q_abs_ADDA_autoX[j]=Pol_X[3,j] 
        Q_ext_ADDA_autoX[j]=Pol_X[1,j]
        
        C_abs_ADDA_autoY[j]=Pol_Y[2,j] 
        C_ext_ADDA_autoY[j]=Pol_Y[0,j]
        C_abs_ADDA_autoX[j]=Pol_X[2,j] 
        C_ext_ADDA_autoX[j]=Pol_X[0,j]
        
        if Csca_option=="yes":        
            Q_sca_ADDA_intY[j]=Pol_Y[5,j] #these output are shown only if -Csca option is selected and it useful to check that everything is correct.
            Q_sca_ADDA_intX[j]=Pol_X[5,j]
        
            C_sca_ADDA_intY[j]=Pol_Y[4,j] 
            C_sca_ADDA_intX[j]=Pol_X[4,j]
        else:
            Q_sca_ADDA_intY[j]=0
            Q_sca_ADDA_intX[j]=0
        
            C_sca_ADDA_intY[j]=0 
            C_sca_ADDA_intX[j]=0
        
    Q_abs_ADDA_auto=(Q_abs_ADDA_autoY+Q_abs_ADDA_autoX)/2
    Q_ext_ADDA_auto=(Q_ext_ADDA_autoY+Q_ext_ADDA_autoX)/2    
    Q_sca_ADDA_int=(Q_sca_ADDA_intY+Q_sca_ADDA_intX)/2
    
    C_abs_ADDA_auto=(C_abs_ADDA_autoY+C_abs_ADDA_autoX)/2
    C_ext_ADDA_auto=(C_ext_ADDA_autoY+C_ext_ADDA_autoX)/2    
    C_sca_ADDA_int=(C_sca_ADDA_intY+C_sca_ADDA_intX)/2
        
        #load Qext directly from CrossSection calculation of output folder, for each wavelength.

        
    Q_ext_ADDA= np.zeros(len(wavelength))
    Q_sca_for_ADDA=np.zeros(len(wavelength))
    Q_sca_back_ADDA=np.zeros(len(wavelength))
    C_ext_ADDA= np.zeros(len(wavelength))    
    C_sca_for_ADDA=np.zeros(len(wavelength))
    C_sca_back_ADDA=np.zeros(len(wavelength))
    Q_sca_Mie_for=np.zeros(len(wavelength))
    Q_sca_Mie_back=np.zeros(len(wavelength))
    Q_ext_Mie=np.zeros(len(wavelength))
    C_sca_ADDA_int_tot=np.zeros(len(wavelength))
    C_sca_ADDA_int_tot_for=np.zeros(len(wavelength))
    C_sca_ADDA_int_tot_back=np.zeros(len(wavelength))
    x=np.zeros(len(wavelength))
    
    for k in range(len(wavelength)):
        #print(wavelength[k])
        x[k]=2*np.pi*a/wavelength[k]    #compute size parameter
        
        dtheta=(S11_theta[len(phi_grid),0]-S11_theta[0,0])*np.pi/180 #step of theta, as defined in muel_scatgrid        
        dphi=(S11_phi[1,0]-S11_phi[0,0])*np.pi/180 #step of phi, as defined in muel_scatgrid
        theta_90_index=int(len(phi_grid)*((len(theta_grid)-1))/2)  #defines the number of row corresponding to theta=90°
        
        if k==0:
            print(f'90° is set as {theta_90_index}')
        
        if substrate_presence=="no": #while working in substrate mode, theta angle must be reversed to assume that positive z axis direction is the impinging light direction...
            C_sca_ADDA_int_tot_for[k] = (((wavelength[k]/1000)**2)/(4*np.pi**2))*np.sum(S11_muel_grid[:theta_90_index,k]*np.sin(S11_theta[:theta_90_index,k]*np.pi/180)*dtheta*dphi)
            C_sca_ADDA_int_tot_back[k] = (((wavelength[k]/1000)**2)/(4*np.pi**2))*np.sum(S11_muel_grid[theta_90_index:,k]*np.sin(S11_theta[theta_90_index:,k]*np.pi/180)*dtheta*dphi)
            C_sca_ADDA_int_tot[k] = (((wavelength[k]/1000)**2)/(4*np.pi**2))*np.sum(S11_muel_grid[:,k]*np.sin(S11_theta[:,k]*np.pi/180)*dtheta*dphi)  #total cross section, assuming no substrate

        if substrate_presence=="yes": #while working in substrate mode, theta angle must be reversed to assume that positive z axis direction is the impinging light direction...
            C_sca_ADDA_int_tot_for[k] = (((wavelength[k]/1000)**2)/(4*np.pi**2))*np.sum(S11_muel_grid[theta_90_index:,k]*np.sin(S11_theta[theta_90_index:,k]*np.pi/180)*dtheta*dphi)/sub_n[k]**2
            C_sca_ADDA_int_tot_back[k] = (((wavelength[k]/1000)**2)/(4*np.pi**2))*np.sum(S11_muel_grid[:theta_90_index,k]*np.sin(S11_theta[:theta_90_index,k]*np.pi/180)*dtheta*dphi)
            C_sca_ADDA_int_tot[k] =  C_sca_ADDA_int_tot_back[k]+C_sca_ADDA_int_tot_for[k]
        
        if size_type=="Size Along X axis":
            Q_ext_Mie[k]=Qext_Mie(PS_n[k]-1j*PS_k[k],wavelength[k],a/2)
            Q_sca_Mie_for[k]=QSca_Mie_for(PS_n[k]-1j*PS_k[k],wavelength[k],a/2)
            Q_sca_Mie_back[k]=QSca_Mie_back(PS_n[k]-1j*PS_k[k],wavelength[k],a/2)
        else:
            Q_ext_Mie[k]=Qext_Mie(PS_n[k]-1j*PS_k[k],wavelength[k],a)
            Q_sca_Mie_for[k]=QSca_Mie_for(PS_n[k]-1j*PS_k[k],wavelength[k],a)
            Q_sca_Mie_back[k]=QSca_Mie_back(PS_n[k]-1j*PS_k[k],wavelength[k],a)    
        
    Q_sca_for_ADDA=C_sca_ADDA_int_tot_for/particle_cross_section
    Q_sca_back_ADDA=C_sca_ADDA_int_tot_back/particle_cross_section
        
    print('Upload from ADDA successfull')
    
    #stacks results and save them inside SimResfolder
    
    resultQ=np.column_stack((wavelength,Q_ext_ADDA_auto,Q_abs_ADDA_auto,Q_ext_ADDA,Q_sca_for_ADDA,Q_sca_back_ADDA,Q_sca_ADDA_int,Q_ext_Mie,Q_sca_Mie_for,Q_sca_Mie_back))
    resultC=np.column_stack((wavelength,C_ext_ADDA_auto,C_abs_ADDA_auto,C_sca_ADDA_int,C_sca_ADDA_int_tot,C_sca_ADDA_int_tot_for,C_sca_ADDA_int_tot_back))    
    
    
    filename_Q = f"Simulation_result_Q_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_fromgrid_{Flag}.txt"
    full_pathQ = os.path.join(SimResfolder, filename_Q)
    np.savetxt(full_pathQ,resultQ,fmt='%.6f',header='wavelength Q_ext_ADDA_auto Q_abs_ADDA_auto Q_ext_ADDA Q_sca_for_ADDA Q_sca_back_ADDA Q_sca_ADDA_int Q_ext_Mie Q_sca_Mie_for Q_sca_Mie_back')
    
    filename_C = f"Simulation_result_C_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_fromgrid_{Flag}.txt"
    full_pathC = os.path.join(SimResfolder, filename_C)
    np.savetxt(full_pathC,resultC,fmt='%.6f',header='wavelength C_ext_ADDA_auto C_abs_ADDA_auto C_sca_ADDA_int C_sca_ADDA_int_tot C_sca_ADDA_int_tot_for C_sca_ADDA_int_tot_back')
  
    print(f'txt files saved in {SimResfolder}')      
    return
    
def plot(a):    
    a_value=str(a)
    
    data_input_C= f"Simulation_result_C_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_{Flag}.txt"
    full_pathC = os.path.join(SimResfolder, data_input_C)
    
    data_input_Q= f"Simulation_result_Q_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_{Flag}.txt"
    full_pathQ = os.path.join(SimResfolder, data_input_Q)
    
    #load efficiencies

    Q_ext_ADDA_auto=np.loadtxt(full_pathQ,skiprows=1, usecols=1)
    Q_abs_ADDA_auto=np.loadtxt(full_pathQ,skiprows=1, usecols=2)
    Q_sca_for_ADDA=np.loadtxt(full_pathQ,skiprows=1, usecols=4)
    Q_sca_back_ADDA=np.loadtxt(full_pathQ,skiprows=1, usecols=5)    
    Q_ext_Mie=np.loadtxt(full_pathQ,skiprows=1, usecols=7)
    Q_sca_Mie_for=np.loadtxt(full_pathQ,skiprows=1, usecols=8)
    Q_sca_Mie_back=np.loadtxt(full_pathQ,skiprows=1, usecols=9)    
    
    Q_ext_ADDA=Q_ext_ADDA_auto  #skip the definition from muel matrx (#TODO: FIX) and use directly the CrossSec output
       
    #load scattering cross sections
    C_ext_ADDA_auto=np.loadtxt(full_pathC,skiprows=1, usecols=1)
    C_abs_ADDA_auto=np.loadtxt(full_pathC,skiprows=1, usecols=2)
    C_sca_ADDA_int=np.loadtxt(full_pathC,skiprows=1, usecols=3)
    C_sca_ADDA_int_tot=np.loadtxt(full_pathC,skiprows=1, usecols=4)
    C_sca_ADDA_int_tot_for=np.loadtxt(full_pathC,skiprows=1, usecols=5)
    C_sca_ADDA_int_tot_back=np.loadtxt(full_pathC,skiprows=1, usecols=6)
    
   
    fig, (ax1)=pylab.subplots(1,1,figsize=(4,8), dpi=600)
    
    ax1.plot(wavelength,Q_ext_Mie,'--',color=(0,1,0),label='ext Mie')
    ax1.plot(wavelength,Q_sca_Mie_for,'--',color=(1,0,0),label='sca for Mie')
    ax1.plot(wavelength,Q_sca_Mie_back,'--',color=(0,0.3,1),label='sca back Mie')
    ax1.plot(wavelength,Q_ext_ADDA,color='green',label='ext ADDA auto')
    ax1.plot(wavelength,Q_sca_for_ADDA,color=(0.5,0,0),label='sca for ADDA')
    ax1.plot(wavelength,Q_sca_back_ADDA,color=(0,0,0.5),label='sca back ADDA')
    ax1.legend(loc="upper right", fontsize="small", frameon=True,framealpha=0)
    ax1.set_title(f"Efficency factors ADDA vs Mie, R={a_value}nm")


    fig2, (ax4)=pylab.subplots(1,1,figsize=(4,8), dpi=600)
    ax4.plot(wavelength,C_sca_ADDA_int_tot_for,color=(0,0,0),label='Csca for')
    ax4.plot(wavelength,C_sca_ADDA_int_tot_back,color=(1,0,0),label='Csca back')
    ax4.plot(wavelength,C_sca_ADDA_int_tot,color=(0,0.8,0),label='Csca for+Csca back ')
    ax4.plot(wavelength,C_sca_ADDA_int,'--',color=(0,0,1),label='Csca out ADDA')
    ax4.plot(wavelength,C_abs_ADDA_auto,color=(1,0,1),label='Cabs auto')
    ax4.plot(wavelength,C_ext_ADDA_auto,color=(0,0.8,1),label='Cext auto')
    ax4.legend(loc="upper right", fontsize="small", frameon=True,framealpha=0)
    ax4.set_title(f"ADDA Scattering cross sections, R={a_value}nm")
      
    pylab.show()
 
def plot_grid(a):    
     a_value=str(a)
     
     data_input_C= f"Simulation_result_C_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_fromgrid_{Flag}.txt"
     full_pathC = os.path.join(SimResfolder, data_input_C)
     
     data_input_Q= f"Simulation_result_Q_R{str(a)}_{int(wavelength[0])}_{int(wavelength[-1])}_{particle_shape}_sub{substrate_presence}_fromgrid_{Flag}.txt"
     full_pathQ = os.path.join(SimResfolder, data_input_Q)
     
     #load efficiencies

     Q_ext_ADDA_auto=np.loadtxt(full_pathQ,skiprows=1, usecols=1)
     Q_abs_ADDA_auto=np.loadtxt(full_pathQ,skiprows=1, usecols=2)
     Q_sca_for_ADDA=np.loadtxt(full_pathQ,skiprows=1, usecols=4)
     Q_sca_back_ADDA=np.loadtxt(full_pathQ,skiprows=1, usecols=5)    
     Q_ext_Mie=np.loadtxt(full_pathQ,skiprows=1, usecols=7)
     Q_sca_Mie_for=np.loadtxt(full_pathQ,skiprows=1, usecols=8)
     Q_sca_Mie_back=np.loadtxt(full_pathQ,skiprows=1, usecols=9)    
     
     Q_ext_ADDA=Q_ext_ADDA_auto  #skip the definition from muel matrx (#TODO: FIX) and use directly the CrossSec output
        
     #load scattering cross sections
     C_ext_ADDA_auto=np.loadtxt(full_pathC,skiprows=1, usecols=1)
     C_abs_ADDA_auto=np.loadtxt(full_pathC,skiprows=1, usecols=2)
     C_sca_ADDA_int=np.loadtxt(full_pathC,skiprows=1, usecols=3)
     C_sca_ADDA_int_tot=np.loadtxt(full_pathC,skiprows=1, usecols=4)
     C_sca_ADDA_int_tot_for=np.loadtxt(full_pathC,skiprows=1, usecols=5)
     C_sca_ADDA_int_tot_back=np.loadtxt(full_pathC,skiprows=1, usecols=6)
     
    
     fig, (ax1)=pylab.subplots(1,1,figsize=(4,6), dpi=600)
     
     ax1.plot(wavelength,Q_ext_Mie,'--',color=(0,1,0),label='ext Mie')
     ax1.plot(wavelength,Q_sca_Mie_for,'--',color=(1,0,0),label='sca for Mie')
     ax1.plot(wavelength,Q_sca_Mie_back,'--',color=(0,0.3,1),label='sca back Mie')
     ax1.plot(wavelength,Q_ext_ADDA,color='green',label='ext ADDA auto')
     ax1.plot(wavelength,Q_sca_for_ADDA,color=(0.5,0,0),label='sca for ADDA')
     ax1.plot(wavelength,Q_sca_back_ADDA,color=(0,0,0.5),label='sca back ADDA')
     ax1.set_xlabel("Wavelength (nm)",fontsize=12, fontweight="bold")
     ax1.set_ylabel("Efficency factor",fontsize=12, fontweight="bold")
     ax1.legend(loc="upper right", fontsize="small", frameon=True,framealpha=0)
     ax1.set_title(f"Efficency factors ADDA vs Mie, R={a_value}nm")


     fig2, (ax4)=pylab.subplots(1,1,figsize=(4,6), dpi=600)
     ax4.plot(wavelength,C_sca_ADDA_int_tot_for,color=(0,0,0),label='Csca for')
     ax4.plot(wavelength,C_sca_ADDA_int_tot_back,color=(1,0,0),label='Csca back')
     ax4.plot(wavelength,C_sca_ADDA_int_tot,color=(0,0.8,0),label='Csca for+Csca back ')
     ax4.plot(wavelength,C_sca_ADDA_int,'--',color=(0,0,1),label='Csca out ADDA')
     ax4.plot(wavelength,C_abs_ADDA_auto,color=(1,0,1),label='Cabs auto')
     ax4.plot(wavelength,C_ext_ADDA_auto,color=(0,0.8,1),label='Cext auto')
     ax4.legend(loc="upper right", fontsize="small", frameon=True,framealpha=0)
     ax4.set_xlabel("Wavelength (nm)",fontsize=12, fontweight="bold")
     ax4.set_ylabel("Cross section (µm²)",fontsize=12, fontweight="bold")
     ax4.set_title(f"ADDA Scattering cross sections, R={a_value}nm")
       
     pylab.show()


def kill_folder(a):  #function to kill the simulation folders and clean the working environment.
    if wavelength.shape[0] > 1:
        initial_run = int(wavelength[0])
        step_lambda = int(wavelength[1] - wavelength[0])
        final_run = initial_run + step_lambda * wavelength.shape[0]
    else:
        initial_run = int(wavelength[0])
        step_lambda = 1   # no step if only one point
        final_run = initial_run+step_lambda
    
    for j, wl in enumerate(wavelength):
        i = int(wl)  # wavelength value, used for folder naming
        run = f"{i:03d}"   
    
        print(j)
            
        folder = (
            Mainfolder
            + "run" + run + "_"
            + particle_shape
            + "_sub" + substrate_presence
            + "_g" + number_of_dipoles_in_a_row
            + "_m" + str(PS_n[j])
            + "_R" + str(a)
            + f"_{Flag}" )
        
        if os.path.exists(folder):
            shutil.rmtree(folder)
            print('folder killed:', folder)
        else:
            print('next folder')
        
    return
 


################################## MIE CALCULATIONS ##############################################

def QSca_Mie_for(m,lambda0,a):  #(refractive index, wavelength, radius) computes the forward scattering efficency. to compare ADDA and Mie results for spheres. needs Miepython module in the win64 directory.
    theta_angle_forM = np.arange(start=0, stop=0.5*np.pi+0.00000001,step=np.pi/100000) #defines set of angle for forward scattering (0-90°) 
    x=2*np.pi*a/lambda0    #compute size parameter
    mu_for=np.cos(theta_angle_forM)  #compute cosine of the selected angle range, serves as input of _mie_S1_S2 function
    norm_int = 6 #Fix Wiscombe normalization, necessary to make results consistent.
    try:
        import miepython
  # try to import inside the function
        S1_for,S2_for=miepython._mie_S1_S2(m, x, mu_for,  norm_int) #compute amplitude function for forward scattering
        Q_sca_for = x**(-2)*np.sum(((np.abs(S1_for)**2 + np.abs(S2_for)**2)*np.sin(theta_angle_forM)*np.pi/100000))
    
    except ModuleNotFoundError:
        print("[WARNING] module miepython not found. Using Q_sca_for=0 instead.")
        Q_sca_for = 0
    except AttributeError:
        print("[WARNING] function _mie_S1_S2 not found in module. Using Q_sca_for=0 instead.")
        Q_sca_for = 0
    except Exception as e:
        print(f"[WARNING] failed to compute Q_sca_for: {e}. Using Q_sca_for=0 instead.")
        Q_sca_for = 0
    

    #print(Q_sca_for)
    #pylab.plot(theta_angle_forM,S1_for,color='red')
    return Q_sca_for

def QSca_Mie_back(m,lambda0,a):  #(refractive index, wavelength, radius) computes the backward scattering efficency. need Miepython module
    theta_angle_forM = np.arange(start=0.5*np.pi, stop=1*np.pi+0.00000001,step=np.pi/100000) #defines set of angle for forward scattering (0-90°) 
    x=2*np.pi*a/lambda0    #compute size parameter
    mu_for=np.cos(theta_angle_forM)  #compute cosine of the selected angle range, serves as input of _mie_S1_S2 function
    norm_int = 6 #Fix Wiscombe normalization, necessary to make results consistent.
    try:
        import miepython  # try to import inside the function
        S1_for,S2_for=miepython._mie_S1_S2(m, x, mu_for,  norm_int) #compute amplitude function for backward scattering
        Q_sca_back = x**(-2)*np.sum(((np.abs(S1_for)**2 + np.abs(S2_for)**2)*np.sin(theta_angle_forM)*np.pi/100000))
    
    except ModuleNotFoundError:
        print("[WARNING] module miepython not found. Using Q_sca_back=0 instead.")
        Q_sca_back = 0
    except AttributeError:
        print("[WARNING] function g_mie_S1_S2 not found in module. Using Q_sca_back=0 instead.")
        Q_sca_back = 0
    except Exception as e:
        print(f"[WARNING] failed to compute Q_sca_back: {e}. Using Q_sca_back=0 instead.")
        Q_sca_back = 0
    
    #print(Q_sca_back)
    #pylab.plot(theta_angle_forM,S1_for,color='red')
    return Q_sca_back

def Qext_Mie(m,lambda0,a): #computes the exctintion scattering efficency. need Miepython module
    x=2*np.pi*a/lambda0
    norm_int=6
    mu_0= np.arange(start=1.0, stop=1.01,step=0.1) #Qext depend only on the value of S for theta=0. mu_0 is an array with one element (costheta=1), because to adjust the normalization the input angle for _mie_s1_s2 must be an array 
    norm_int = 6
    try:
        import miepython  # try to import inside the function
        S1_0, S2_0 = miepython._mie_S1_S2(m,x,mu_0, norm_int) #compute amplitude function at theta=0
        Q_ext_A = 4*(1/x**2)*S1_0.real #compute Qext
    
    except ModuleNotFoundError:
        print("[WARNING] module miepython not found. Using Q_ext=0 instead.")
        Q_ext_A = 0
    except AttributeError:
        print("[WARNING] function _mie_S1_S2 not found in module. Using Q_ext=0 instead.")
        Q_ext_A = 0
    except Exception as e:
        print(f"[WARNING] failed to compute Q_ext_A: {e}. Using Q_ext=0 instead.")
        Q_ext_A = 0
    
    Q_ext=float(Q_ext_A[0])
    return Q_ext


################################# GUI ####################################################

# === global vars that your functions depend on ===
particle_shape = ""
substrate_presence = ""
grid_value = 0
number_of_dipoles_in_a_row = ""
theta_grid = None
phi_grid = None
particle_cross_section = None
parallel_processes = 1

# # === GUI callbacks ===

def run_upload_plot_sim():
    global keep_bat_files,size_type,phi_grid,particle_cross_section,store_scat_grid_option,Csca_option,particle_shape, substrate_presence, grid_value, number_of_dipoles_in_a_row, theta_grid, parallel_processes,ntheta_value,custom_shape,distance_from_surf,Flag
    try:
        grid_value = int(entry_grid.get())
        number_of_dipoles_in_a_row = str(grid_value)
        substrate_presence = substrate_var.get()
        a = float(entry_a.get())
        particle_shape = entry_shape.get()
        parallel_processes = int(entry_parallel.get())
        theta_step = float(entry_theta_step.get())
        theta_stop = float(entry_theta_stop.get())
        theta_grid = np.arange(start=0, stop=theta_stop+theta_step/2, step=theta_step)
        ntheta_value=len(theta_grid)-1    
        phi_step = float(entry_phi_step.get())
        phi_stop = float(entry_phi_stop.get())        
        phi_grid = np.arange(start=0, stop=phi_stop+phi_step/2, step=phi_step)
        particle_cross_section = np.pi * (a/1000)**2
        
        custom_shape = custom_shape_var.get()
        distance_from_surf=distance_from_surf_var.get()
        Flag=entry_Flag.get()
        Csca_option=Csca_option_var.get()
        store_scat_grid_option=store_scat_grid_option_var.get()
        size_type=size_type_var.get()
        keep_bat_files=keep_bat_files_var.get()
        
        if store_scat_grid_option=="yes":
            print('Performing simulation on whole Theta-Phi grid')
            go(a)
            upload_grid(a)
            plot_grid(a)
        if store_scat_grid_option=="no":
            print('Performing simulation on just the Theta grid')
            go(a)
            upload(a)
            plot(a)
            
        messagebox.showinfo("Done", "Theta-only simulation finished.")
    except Exception as e:
        messagebox.showerror("Error", str(e))

def run_sim():
    global keep_bat_files,size_type,phi_grid,particle_cross_section,store_scat_grid_option,Csca_option,particle_shape, substrate_presence, grid_value, number_of_dipoles_in_a_row, theta_grid, parallel_processes,ntheta_value,custom_shape,distance_from_surf,Flag
    try:
        grid_value = int(entry_grid.get())
        number_of_dipoles_in_a_row = str(grid_value)
        substrate_presence = substrate_var.get()
        a = float(entry_a.get())
        particle_shape = entry_shape.get()
        parallel_processes = int(entry_parallel.get())
        theta_step = float(entry_theta_step.get())
        theta_stop = float(entry_theta_stop.get())
        theta_grid = np.arange(start=0, stop=theta_stop+theta_step/2, step=theta_step)
        ntheta_value=len(theta_grid)-1    
        phi_step = float(entry_phi_step.get())
        phi_stop = float(entry_phi_stop.get())        
        phi_grid = np.arange(start=0, stop=phi_stop+phi_step/2, step=phi_step)
        particle_cross_section = np.pi * (a/1000)**2
        
        custom_shape = custom_shape_var.get()
        distance_from_surf=distance_from_surf_var.get()
        Flag=entry_Flag.get()
        Csca_option=Csca_option_var.get()
        store_scat_grid_option=store_scat_grid_option_var.get()
        size_type=size_type_var.get()
        keep_bat_files=keep_bat_files_var.get()
        
        if store_scat_grid_option=="yes":
            print('Performing simulation on whole Theta-Phi grid')
            go(a)
        if store_scat_grid_option=="no":
            print('Performing simulation on just the Theta grid')
            go(a)
        
        messagebox.showinfo("Done", "Theta-only simulation finished.")
    except Exception as e:
        messagebox.showerror("Error", str(e))
    
def upload_sim():
    global size_type,phi_grid,particle_cross_section,store_scat_grid_option,Csca_option,particle_shape, substrate_presence, grid_value, number_of_dipoles_in_a_row, theta_grid, parallel_processes,ntheta_value,custom_shape,distance_from_surf,Flag
    try:
        grid_value = int(entry_grid.get())
        number_of_dipoles_in_a_row = str(grid_value)
        substrate_presence = substrate_var.get()
        a = float(entry_a.get())
        particle_shape = entry_shape.get()
        parallel_processes = int(entry_parallel.get())
        theta_step = float(entry_theta_step.get())
        theta_stop = float(entry_theta_stop.get())
        theta_grid = np.arange(start=0, stop=theta_stop+theta_step/2, step=theta_step)
        ntheta_value=len(theta_grid)-1    
        phi_step = float(entry_phi_step.get())
        phi_stop = float(entry_phi_stop.get())        
        phi_grid = np.arange(start=0, stop=phi_stop+phi_step/2, step=phi_step)
        particle_cross_section = np.pi * (a/1000)**2
        size_type=size_type_var.get()
        
        custom_shape = custom_shape_var.get()
        distance_from_surf=distance_from_surf_var.get()
        Flag=entry_Flag.get()
        Csca_option=Csca_option_var.get()
        store_scat_grid_option=store_scat_grid_option_var.get()
        
        if store_scat_grid_option=="yes":
            upload_grid(a)
        if store_scat_grid_option=="no":
            upload(a)
        
        messagebox.showinfo("Done", "Upload finished.")
    except Exception as e:
        messagebox.showerror("Error", str(e))


def plot_sim():
    global phi_grid,particle_cross_section,store_scat_grid_option,Csca_option,particle_shape, substrate_presence, grid_value, number_of_dipoles_in_a_row, theta_grid, parallel_processes,ntheta_value,custom_shape,distance_from_surf,Flag
    try:
        grid_value = int(entry_grid.get())
        number_of_dipoles_in_a_row = str(grid_value)
        substrate_presence = substrate_var.get()
        a = float(entry_a.get())
        particle_shape = entry_shape.get()
        parallel_processes = int(entry_parallel.get())
        theta_step = float(entry_theta_step.get())
        theta_stop = float(entry_theta_stop.get())
        theta_grid = np.arange(start=0, stop=theta_stop+theta_step/2, step=theta_step)
        ntheta_value=len(theta_grid)-1    
        phi_step = float(entry_phi_step.get())
        phi_stop = float(entry_phi_stop.get())        
        phi_grid = np.arange(start=0, stop=phi_stop+phi_step/2, step=phi_step)
        particle_cross_section = np.pi * (a/1000)**2
        
        
        custom_shape = custom_shape_var.get()
        distance_from_surf=distance_from_surf_var.get()
        Flag=entry_Flag.get()
        Csca_option=Csca_option_var.get()
        store_scat_grid_option=store_scat_grid_option_var.get()
        ntheta_value=len(theta_grid)-1    

        
        if store_scat_grid_option=="yes":
            plot_grid(a)
        if store_scat_grid_option=="no":
            plot(a)

        
        #messagebox.showinfo("Done", "Plot done.")
    except Exception as e:
        messagebox.showerror("Error", str(e))


def kill_folder_gui():
    global Flag,particle_shape, substrate_presence, grid_value, number_of_dipoles_in_a_row, theta_grid, phi_grid, particle_cross_section, parallel_processes
    try:
        grid_value = int(entry_grid.get())
        number_of_dipoles_in_a_row = str(grid_value)
        substrate_presence = substrate_var.get()
        a = float(entry_a.get())
        particle_shape = entry_shape.get()
        parallel_processes = int(entry_parallel.get())
        theta_step = float(entry_theta_step.get())
        theta_stop = float(entry_theta_stop.get())
        phi_step = float(entry_phi_step.get())
        phi_stop = float(entry_phi_stop.get())
        Flag=entry_Flag.get()
        
        theta_grid = np.arange(start=0, stop=theta_stop+theta_step/2, step=theta_step)
        phi_grid = np.arange(start=0, stop=phi_stop+phi_step/2, step=phi_step)
        particle_cross_section = np.pi * (a/1000)**2

        kill_folder(a)
        messagebox.showinfo("Done", "Folders killed.")
    except Exception as e:
        messagebox.showerror("Error", str(e))


# === build GUI ===
root = tk.Tk()
root.title("ADDA Simulation GUI")

# parallel processes
tk.Label(root, text="Parallel Processes:").grid(row=0, column=0, sticky="w")
entry_parallel = tk.Entry(root, width=6)
entry_parallel.insert(0, "4")
entry_parallel.grid(row=0, column=1)

# grid value
tk.Label(root, text="Grid Value:").grid(row=1, column=0, sticky="w")
entry_grid = tk.Entry(root, width=6)
entry_grid.insert(0, "32")
entry_grid.grid(row=1, column=1)

# particle size
tk.Label(root, text="Particle dimension [nm]:").grid(row=2, column=0, sticky="w")
entry_a = tk.Entry(root, width=6)
entry_a.insert(0, "50")
entry_a.grid(row=2, column=1)

# Define the particle dimension as equivalent disc radius or size along x axis.
tk.Label(root, text="Particle dimension type:").grid(row=3, column=0, sticky="w")
size_type_var = tk.StringVar(value="Equivalent radius")
ttk.Combobox(root, textvariable=size_type_var, values=["Equivalent radius", "Size Along X axis"], width=20).grid(row=3, column=1)

# substrate presence
tk.Label(root, text="Substrate:").grid(row=4, column=0, sticky="w")
substrate_var = tk.StringVar(value="no")
ttk.Combobox(root, textvariable=substrate_var, values=["yes", "no"], width=5).grid(row=4, column=1)

# distance from surface 
tk.Label(root, text="Distance between particle center and surface [nm]:").grid(row=5, column=0, sticky="w")
distance_from_surf_var = tk.Entry(root, width=6)
distance_from_surf_var.insert(0, "50")
distance_from_surf_var.grid(row=5, column=1)

# --- Custom Shape (Yes/No) ---
tk.Label(root, text="Custom Shape: ").grid(row=6, column=0, sticky="w", pady=2)
custom_shape_var = tk.StringVar(value="no")  # default value
ttk.Combobox(root, textvariable=custom_shape_var, values=["yes", "no"], width=5).grid(row=6, column=1)

# particle shape
tk.Label(root, text="Particle Shape:").grid(row=7, column=0, sticky="w")
entry_shape = tk.Entry(root, width=20)
entry_shape.insert(0, "sphere")
entry_shape.grid(row=7, column=1)

# theta grid
tk.Label(root, text="Theta Stop [deg]:").grid(row=8, column=0, sticky="w")
entry_theta_stop = tk.Entry(root, width=6)
entry_theta_stop.insert(0, "180")
entry_theta_stop.grid(row=8, column=1)

tk.Label(root, text="Theta Step [deg]:").grid(row=9, column=0, sticky="w")
entry_theta_step = tk.Entry(root, width=6)
entry_theta_step.insert(0, "0.5")
entry_theta_step.grid(row=9, column=1)

# phi grid
tk.Label(root, text="Phi Stop [deg]:").grid(row=10, column=0, sticky="w")
entry_phi_stop = tk.Entry(root, width=6)
entry_phi_stop.insert(0, "360")
entry_phi_stop.grid(row=10, column=1)

tk.Label(root, text="Phi Step [deg]:").grid(row=11, column=0, sticky="w")
entry_phi_step = tk.Entry(root, width=6)
entry_phi_step.insert(0, "2")
entry_phi_step.grid(row=11, column=1)

tk.Label(root, text="Flag:").grid(row=12, column=0, sticky="w")
entry_Flag = tk.Entry(root, width=20)
entry_Flag.insert(0, "")
entry_Flag.grid(row=12, column=1)

tk.Label(root, text="Internal Csca integration:").grid(row=13, column=0, sticky="w")
Csca_option_var = tk.StringVar(value="yes")
ttk.Combobox(root, textvariable=Csca_option_var, values=["yes", "no"], width=5).grid(row=13, column=1)

tk.Label(root, text="Enable Phi scattering grid:").grid(row=14, column=0, sticky="w")
store_scat_grid_option_var = tk.StringVar(value="yes")
ttk.Combobox(root, textvariable=store_scat_grid_option_var, values=["yes", "no"], width=5).grid(row=14, column=1)



tk.Label(root, text="Keep .bat files (debug option):").grid(row=15, column=0, sticky="w")
keep_bat_files_var = tk.StringVar(value="no")
ttk.Combobox(root, textvariable=keep_bat_files_var, values=["yes", "no"], width=5).grid(row=15, column=1)

# buttons

tk.Button(root, text="Run,upload and plot", command=run_upload_plot_sim).grid(row=16, column=0, columnspan=2, padx=10, pady=5)
tk.Button(root, text="Run", command=run_sim).grid(row=17, column=0, columnspan=2, padx=10, pady=5)
tk.Button(root, text="Upload", command=upload_sim).grid(row=18, column=0, columnspan=2, padx=10, pady=5)
tk.Button(root, text="Plot", command=plot_sim).grid(row=19, column=0, columnspan=2, padx=10, pady=5)

tk.Button(root, text="Kill simulation folders", command=kill_folder_gui).grid(row=20, column=0, columnspan=2, padx=10, pady=10)

root.mainloop()






