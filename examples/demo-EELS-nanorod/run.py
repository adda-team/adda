import sys, os, multiprocessing
if sys.path[1] != os.path.abspath(__file__ + "/../../") : sys.path.insert(1,os.path.abspath(__file__ + "/../../")) #This is to import ADDA Wrapper from parent directory
import ADDAWrapper as aw

# PythonADDAWrapper parameters
aw_parameters = dict(
    #adda_exec = "../../win64/adda.exe", #path to ADDA executable
    adda_exec = os.path.abspath(__file__ + "/../../../src/seq/adda"), #path to ADDA executable
    parallel_procs = multiprocessing.cpu_count()-1, #number of parallel processes is equal to the number of processor cores minus 1

    mp_file = os.path.abspath(__file__ + "/../../refractive_index/" + "Au_JHW.csv"), #file with refractive index of the particle, each string contains: ev,mp_re,mp_im
    ev_range = (0.5,3), #[eV]. Used in "spectrum_" functions. (ev_min,ev_max): range from ev_min[eV] to ev_max[eV]
    ev = 1.95, #[eV]. Used in "scan_" and "extrapolation_" functions.
    
    spectrumline_startpoint = (10,0), # (x,y) [nm]
    spectrumline_endpoint = (10,80), # (x,y) [nm]
    spectrumline_points = 15, #how many points, including startpoint and endpoint
    
    #Used in "scan_" functions. Beam propagation must be orthogonal to the grid.
    #So "prop" must be "0 0 whatever" and rotations with "orient" must be made by 90 degrees.
    scan_x_range = (-15,15), #[nm], (x_left, x_right)
    scan_y_range = (-50,50), #[nm], (y_bottom, y_top)
    scan_step = 1 #dipoles per each step, must be an integer >= 1. Set to 1 for the finest resolution.
    #The beam must always blast exactly in the middle between the dipoles,
    #so start and stop coordinates will be adjusted, covering more area than you entered. Obligatory to use "no_vol_cor" with scan.
)

# Not an arg yet, soon to be implemented into ADDA
hd = (92.6-7.8)/7.8

# ADDA command line arguments
adda_cmdlineargs = dict(
    # Particle
    shape = f"capsule {hd}",
    size = 7.8, #[nm]
    grid = 8, #dipoles per x-axis size of the particle
    mhost = "1.45 0", #refractive index of the host medium
    orient = "90 90 0", #rotating the particle
    
    # Beam
    beam = "electron 100 10 20 0", #Energy[keV] x[nm] y[nm] z[nm] m_host
    prop = "0 0 -1", #beam propagation direction vector
    
    # Precision and performance
    eps = 4, #Residual norm
    
    # Additional options
    sym = "enf", #Do not simulate second polarization
    scat_matr = "none", #Do not calculate the Mueller matrix
    no_vol_cor = "", #Disable volume correction
    iter = "qmr2", #Iterative solver
    pol = "igt_so", #Polarizability prescription
    int = "igt 3", #Interaction term
)

### Executing commands
if __name__ == '__main__': 

    # Execute spectra simulations for different positions of the beam to find plasmon peaks
    dirname = os.path.abspath(__file__ + "/../" + "spectrumline")
    aw.spectrumline_execute(aw_parameters,adda_cmdlineargs,dirname)
    # Collect and plot EELS spectra
    aw.spectrumline_collect("Peels",dirname)
    aw.spectrumline_plot("Peels",dirname)
    aw.spectrumline_collect("Pcl",dirname)
    aw.spectrumline_plot("Pcl",dirname)
    
    # Execute scan of particle's cross-section for single energy ev
    dirname = dirname = os.path.abspath(__file__ + "/../" + "scan")
    aw.scan_execute(aw_parameters,adda_cmdlineargs,dirname)
    # Collect and map scanned EELS probabilities on particle's cross-section
    aw.scan_collect("Peels",dirname)
    aw.scan_plot("Peels",dirname)
    aw.scan_collect("Pcl",dirname)
    aw.scan_plot("Pcl",dirname)


