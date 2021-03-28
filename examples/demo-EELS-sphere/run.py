import sys, os, multiprocessing
if sys.path[1] != os.path.abspath(__file__ + "/../../") : sys.path.insert(1,os.path.abspath(__file__ + "/../../")) #This is to import ADDA Wrapper from parent directory
import ADDAWrapper as aw

# PythonADDAWrapper parameters
aw_parameters = dict(
    #adda_exec = "../../win64/adda.exe", #path to ADDA executable
    adda_exec = os.path.abspath(__file__ + "/../../../src/seq/adda"), #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-1, #number of parallel processes is equal to the number of processor cores minus 1

    mp_file = os.path.abspath(__file__ + "/../../refractive_index/" + "Ag_JHW.csv"), #file with refractive index of the particle, each string contains: ev,mp_re,mp_im
    ev_range = (2.5,4.5), #[eV]. Used in "spectrum_" functions. (ev_min,ev_max): range from ev_min[eV] to ev_max[eV]
    ev = 3.45, #[eV]. Used in "scan_" and "extrapolation_" functions.
    
    #Used in "scan_" functions. Beam propagation must be orthogonal to the grid.
    #So "prop" must be "0 0 whatever" and rotations with "orient" must be made by 90 degrees.
    scan_x_range = (0,30), #[nm], (x_left, x_right)
    scan_y_range = (0,30), #[nm], (y_bottom, y_top)
    scan_step = 1 #dipoles per each step, must be an integer >= 1
    #The beam must always blast exactly in the middle between the dipoles,
    #so start and stop coordinates will be adjusted, covering more area than you entered. Obligatory to use "no_vol_cor" with scan.
)

# ADDA command line arguments
adda_cmdlineargs = dict(
    # Particle
    shape = "sphere",
    size = 40, #[nm]
    grid = 32, #dipoles per x-axis size of the particle
    mhost = "1 0", #refractive index of the host medium
    
    # Beam
    beam = "electron 100 60 0 0", #Energy[keV] x[nm] y[nm] z[nm] m_host
    prop = "0 0 -1", #beam propagation direction vector
    
    # Precision and performance
    eps = 4, #Residual norm
    
    # Additional options
    sym = "enf", #Do not simulate second polarization
    scat_matr = "none", #Do not calculate the Mueller matrix
    no_vol_cor = "", #Disable volume correction
    iter = "qmr2", #Iterative solver
    pol = "igt_so", #Polarizability prescription
    int = "igt 5", #Interaction term
)

### Executing commands
if __name__ == '__main__': 
    
    # Execute spectrum
    dirname = os.path.abspath(__file__ + "/../" + "spectrum")
    aw.spectrum_execute(aw_parameters,adda_cmdlineargs,dirname)
    # Collect and plot EELS probabilities
    aw.spectrum_collect("Peels",dirname)
    aw.spectrum_plot("Peels",dirname)
    # Collect and plot CL probabilities
    aw.spectrum_collect("Pcl",dirname)
    aw.spectrum_plot("Pcl",dirname)
    
    # Execute extrapolation for single energy ev
    dirname = os.path.abspath(__file__ + "/../" + "extrapolation")
    aw.extrapolation_execute(aw_parameters,adda_cmdlineargs,dirname)
    # Collect and plot EELS probabilities for different y~1/grid with extrapolated value in y=0 (with errorbar)
    aw.extrapolation_collect("Peels",dirname)
    aw.extrapolation_plot("Peels",dirname)
    
    # Execute spectrum with extrapolation at each energy ev
    dirname = os.path.abspath(__file__ + "/../" + "spectrum_with_extrapolation")
    aw.spectrum_with_extrapolation_execute(aw_parameters,adda_cmdlineargs,dirname)
    # Collect and plot fit EELS probabilities spectrum (with errorbar)
    aw.spectrum_with_extrapolation_collect("Peels",dirname)
    aw.spectrum_with_extrapolation_plot("Peels",dirname)
    
    # Execute scan of particle's cross-section for single energy ev
    dirname = dirname = os.path.abspath(__file__ + "/../" + "scan")
    aw.scan_execute(aw_parameters,adda_cmdlineargs,dirname)
    # Collect and map scanned EELS probabilities on particle's cross-section
    aw.scan_collect("Peels",dirname)
    aw.scan_plot("Peels",dirname)
    aw.scan_collect("Pcl",dirname)
    aw.scan_plot("Pcl",dirname)


