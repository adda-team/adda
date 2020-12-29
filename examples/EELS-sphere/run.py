import multiprocessing, sys, time
if sys.path[1] != r"../" : sys.path.insert(1,r"../") #This is to import ADDA Wrapper from parent directory
import ADDAWrapper as aw

# PythonADDAWrapper parameters
aw_parameters = dict(
    #adda_exec = "../../win64/adda.exe", #path to ADDA executable
    adda_exec = "../../src/seq/adda", #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-1, #number of parallel processes is equal to the number of processor cores minus 1

    mp_file = "../refractive_index/Ag_JC.csv", #file with refractive index of the particle, each string contains: ev,mp_re,mp_im
    ev_range = (3,5), #[eV]. Used in "spectrum_" functions. (ev_min,ev_max): range from ev_min[eV] to ev_max[eV]
    ev = 3.6, #[eV]. Used in "scan_" and "extrapolation_" functions.
    
    #Used in "scan_" functions. Assuming beam propagation vector = (0,0,whatever).
    scan_x_range = (0,100), #[nm], (x_left, x_right)
    scan_y_range = (0,100), #[nm], (y_bottom, y_top)
    scan_step = 1 #dipoles per each step, must be an integer >= 1
    #The beam must always blast exactly in the middle between the dipoles,
    #so start and stop coordinates will be adjusted, covering more area than you entered. Obligatory to use "no_vol_cor" with scan.
)

# Not an arg yet, soon to be implemented into ADDA
mh = 1 #refractive index of the host medium

# ADDA command line arguments
adda_cmdlineargs = dict(
    # Particle
    shape = "sphere",
    size = 150, #[nm]
    grid = 16, #dipoles per axis
    mh = mh, #refractive index of the host medium
    
    # Beam
    beam = f"electron 100 100 0 0 {mh}", #Energy[keV] x[nm] y[nm] z[nm] m_host
    prop = "0 0 -1", #beam propagation direction vector
    
    # Precision and performance
    eps = 5, #Residual norm
    
    # Additional options
    sym = "enf", #Do not simulate second polarization
    scat_matr = "none", #Do not calculate the Mueller matrix
    pol = "cm", #Polarizability prescription
    no_vol_cor = "" #Disable volume correction
)

### Executing commands
if __name__ == '__main__': 
    
    # Execute spectrum
    dirname="spectrum"
    aw.spectrum_execute(aw_parameters,adda_cmdlineargs,dirname)
    # Collect and plot EELS probabilities
    aw.spectrum_collect("Peels",dirname)
    aw.spectrum_plot("Peels",dirname)
    # Collect and plot CL probabilities
    aw.spectrum_collect("Pcl",dirname)
    aw.spectrum_plot("Pcl",dirname)
    
    # Execute extrapolation for single energy ev
    dirname="extrapolation"
    aw.extrapolation_execute(aw_parameters,adda_cmdlineargs,dirname)
    # Collect and plot EELS probabilities for different y~1/grid with extrapolated value in y=0 (with errorbar)
    aw.extrapolation_collect("Peels",dirname)
    aw.extrapolation_plot("Peels",dirname)
    
    # Execute spectrum with extrapolation at each energy ev
    dirname = "spectrum_with_extrapolation"
    aw.spectrum_with_extrapolation_execute(aw_parameters,adda_cmdlineargs,dirname)
    # Collect and plot fit EELS probabilities spectrum (with errorbar)
    aw.spectrum_with_extrapolation_collect("Peels",dirname)
    aw.spectrum_with_extrapolation_plot("Peels",dirname)
    
    # Execute scan of particle's cross-section for single energy ev
    dirname = "scan"
    aw.scan_execute(aw_parameters,adda_cmdlineargs,dirname)
    # Collect and map scanned EELS probabilities on particle's cross-section
    aw.scan_collect("Peels",dirname)
    aw.scan_plot("Peels",dirname)


