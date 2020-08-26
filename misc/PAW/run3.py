'''
noname v0.3
Script performs batch simulations using ADDA. For more info read the following articles.
ADDA: Maxim A. Yurkin, "User Manual for the Discrete Dipole Approximation Code ADDA"
EELS and CL: Kichigin&Yurkin, "Simulating Electron-energy-loss Spectroscopy and Cathodoluminescence with the Discrete Dipole Approximation"
Extrapolation: Yurkin et al., "Convergence of the discrete dipole approximation. II. An extrapolation technique to increase the accuracy"
'''
import os, shutil, re, csv, time, multiprocessing, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

cmdline = "/Users/user/Documents/GitHub/adda/src/seq/adda" #path to adda executable
spectrum_path = "spectrum" #default path for storing spectrum data
extrapolation_path = "extrapolation" #default path for storing extrapolation data
spectrum_with_extrapolation_path = "spectrum_with_extrapolation" #default path for storing spectrum_with_extrapolation data
scan_path = "scan" #default path for storing scan data

#Particle parameters
size = 75 #nm
shape = "sphere"
grid = 10 #dipoles per axis
mp_file = "Ag_JC_Garcia1545.csv" #m_particle, each string contains: ev,n,k
mh = 1 #m_host, refractive index of the host medium

#Electron beam parameters
e_energy = 100 #keV
r0 = (0,100,0) #nm, beam position in space
prop = (0,0,-1) #beam propagation direction vector
beam = "electron " + str(e_energy) + " " + str(r0[0]) + " " + str(r0[1]) + " " + str(r0[2]) + " " + str(mh)

#Scan parameters. Assuming prop = (0,0,whatever).
(x_left,x_right) = (0,101*1.1) #nm
(y_bottom,y_top) = (0,20*1.1) #nm
#The beam must always blast exactly in the middle between the dipoles,
#so start and stop coordinates will be adjusted, covering more area than you entered

#Precision and performance
eps = 4 #Residual norm
procsnumber = multiprocessing.cpu_count() #number of parallel processes is equal to the number of processor cores
procsnumber = 3 #manually select number of processes

#Constructing command line with some additional options
cmdline += " -shape " + shape
cmdline += " -size " + str(size)
cmdline += " -eps " + str(eps)
cmdline += " -sym enf" #Do not simulate second polarization
cmdline += " -scat_matr none" #Do not calculate the Mueller matrix
cmdline += " -no_vol_cor" #Disable volume correction
cmdline += " -pol cm" #Polarizability prescription
# cmdline += " -surf 5 2 0" #Surface mode
# cmdline += " -orient 0 90 0" #Rotate particle
# cmdline += " -store_int_field" #Save internal E-field into a file

def exec_cmdline(ev,mre,mim,beam_i,runpath_i):
    lam = 1239.8419842361123824 / (ev * mh)
    cmdline_i = cmdline
    cmdline_i += " -lambda " + str(lam)
    cmdline_i += " -m " + str(mre) + " " + str(mim)
    cmdline_i += " -grid " + str(grid)
    cmdline_i += " -beam " + beam_i
    cmdline_i += " -dir '" + runpath_i + "'"
    cmdline_i += " > /dev/null"
    # print(cmdline_i)
    flag = os.system(cmdline_i)
    if flag != 0:
        print("'" + cmdline_i + "' ran with exit code ", flag)
        return flag
    return 0

def spectrum_execute():
    start_time = time.time()
    shutil.rmtree(spectrum_path, ignore_errors=True)
    os.makedirs(spectrum_path, exist_ok=True)
    print_log("--- Spectrum: executing simulations", path=spectrum_path)
    mdata = np.genfromtxt(mp_file,delimiter=',')
    print_log("mp_file: " + mp_file, path=spectrum_path)
    print_log("From " + str(mdata[0][0]) + " to " + str(mdata[-1][0]) + " eV", path=spectrum_path)
    cmdlineargs = []
    for i in mdata:
        runpath_i = spectrum_path + "/" + str(i[0])
        cmdlineargs.append((i[0],i[1]/mh,i[2]/mh,beam,runpath_i))
    multiprocessing.Pool(procsnumber).starmap(exec_cmdline, cmdlineargs, 1)
    print_log("--- %s seconds" % round((time.time() - start_time),2), path=spectrum_path)

def spectrum_collect(match):
    print("--- Spectrum: collecting results for", match)
    start_time = time.time()
    evs = sorted([float(d.name) for d in os.scandir(spectrum_path) if d.is_dir()])
    values = []
    for ev in evs:
        with open(spectrum_path + "/" + str(ev) + "/CrossSec-Y", "r") as file:
            for line in file:
                if match in line:
                    value = float(re.findall("[ \t][-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?",line)[0])
                    values.append(value)
                    break
    with open(spectrum_path + "/" + match+".csv", 'w') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(["ev",match])
            writer.writerows(zip(evs,values))
            print("Saved to " + spectrum_path + "/"+match+".csv")
    print("--- %s seconds" % round((time.time() - start_time),2))
    
def spectrum_plot(match):
    print("--- Spectrum: plotting results for", match)
    data = np.genfromtxt(spectrum_path + "/"+match+".csv",delimiter=',')[1:]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(data[:,0], data[:,1], label=label_for_plot(match))
    ax.set_xlabel("eV")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(which="both", axis="x", linestyle="dotted")
    ax.legend() 
    fig.savefig(spectrum_path + "/"+match+".png", bbox_inches='tight', dpi=300)

def extrapolation_execute(ev):
    start_time = time.time()
    shutil.rmtree(extrapolation_path, ignore_errors=True)
    os.makedirs(extrapolation_path, exist_ok=True)
    print_log("--- Extrapolation: executing simulations", path=extrapolation_path)
    lam = 1239.8419842361123824 / (ev * mh)
    mdata = np.genfromtxt(mp_file,delimiter=',')
    flag = 0
    for line in mdata:
        if float(line[0]) == float(ev):
            flag = 1
            break
    if flag == 0:
        print("ERROR: ev =",ev,"not found in",mp_file)
        return
    print_log("ev = " + str(ev), path=extrapolation_path)
    print_log("mp_re = " + str(line[1]), path=extrapolation_path)
    print_log("mp_im = " + str(line[2]), path=extrapolation_path)
    mre = float(line[1])/mh
    mim = float(line[2])/mh
    m_abs = math.sqrt(mre**2 + mim**2)
    y_min = (2*math.pi/lam)*(size/grid)*m_abs #y = k*d*|m|
    y_max = 4*y_min
    ys = np.exp(np.linspace(math.log(y_min), math.log(y_max), 9))
    grids = np.rint((2*math.pi/lam)*(size/ys)*m_abs).astype(np.int32)
    print_log("grids: " + str(grids), path=extrapolation_path)
    cmdlineargs = []
    for grid_i in grids:
        runpath_i = extrapolation_path + "/" + str(grid_i)
        cmdlineargs.append((ev,mre,mim,beam,runpath_i))
    multiprocessing.Pool(procsnumber).starmap(exec_cmdline, cmdlineargs, 1)
    with open(extrapolation_path + "/ev.csv", 'w') as file:
        csv.writer(file, delimiter=',').writerow([ev])
    print_log("--- %s seconds" % round((time.time() - start_time),2), path=extrapolation_path)
    
def extrapolation_collect(match, silent=False):
    if silent == False:
        print("--- Extrapolation: collecting results for", match)
    start_time = time.time()
    ev = np.genfromtxt(extrapolation_path + "/ev.csv",delimiter=',')
    lam = 1239.8419842361123824 / (ev * mh)
    mdata = np.genfromtxt(mp_file,delimiter=',')
    flag = 0
    for line in mdata:
        if float(line[0]) == float(ev):
            flag = 1
            break
    if flag == 0:
        print("ERROR: ev =",ev,"not found in",mp_file)
        return
    mre = float(line[1])/mh
    mim = float(line[2])/mh
    m_abs = math.sqrt(mre**2 + mim**2)
    grids = np.array(sorted([int(d.name) for d in os.scandir(extrapolation_path) if d.is_dir()]))
    values = []
    for grid_i in grids:
        with open(extrapolation_path + "/" + str(grid_i) + "/CrossSec-Y", "r") as file:
            for line in file:
                if match in line:
                    value = float(re.findall("[ \t][-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?",line)[0])
                    values.append(value)
                    break
    ys = (2*math.pi/lam)*(size/grids)*m_abs #y = k*d*|m|
    weights = ys**-3
    fit,cov = np.polyfit(ys, values, 2, w=weights, cov=True)
    a = np.flip(fit)
    error = 2*np.sqrt(np.flip(np.diag(cov)))
    with open(extrapolation_path + "/"+match+".csv", 'w') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(["grids","ys","values"])
            writer.writerows(zip(grids,ys,values))
            if silent == False:
                print("Saved to " + extrapolation_path + "/"+match+".csv")
                print(match,"= %.8e +/- %.8e" % (a[0],error[0]))
    with open(extrapolation_path + "/"+match+"_fit.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(["a[i]","error[i]"])
        writer.writerows(zip(a,error))
        if silent == False:
            print("Saved to " + extrapolation_path + "/"+match+"_fit.csv")
    if silent == False:
        print("--- %s seconds" % round((time.time() - start_time),2))
    
def extrapolation_plot(match):
    #Add exact Mie solution
    # miedata = np.genfromtxt("Peels_mie1545.csv",delimiter=',')
    # ev = np.genfromtxt(extrapolation_path + "/ev.csv",delimiter=',')
    # flag = 0
    # for line in miedata:
    #     if float(line[0]) == float(ev):
    #         flag = 1
    #         break
    # if flag == 0:
    #     print("ERROR: ev =",ev,"not found in mie_data")
    #     return
    # plt.plot(0, line[1], label="Peels_Mie", marker="o", color="red")
    
    print("--- Extrapolation: plotting results for", match)
    data = np.genfromtxt(extrapolation_path + "/"+match+".csv",delimiter=',')[1:]
    plt.plot(data[:,1], data[:,2], label=match, marker="o", linestyle="none")
    ys_fitted = np.linspace(data[:,1][0],0,100)
    results_fit = np.genfromtxt(extrapolation_path + "/"+match+"_fit.csv",delimiter=',')[1:]
    a = results_fit[:,0]
    error = results_fit[:,1]
    points_fitted = a[0] + a[1]*ys_fitted + a[2]*ys_fitted**2
    plt.figure()
    plt.plot(ys_fitted, points_fitted, label=match+"_fit", color="black")
    plt.errorbar(0, a[0], yerr=error[0], color="black", marker="s", capsize=3, barsabove=True)
    plt.xlabel("y = kd|m|")
    
    plt.show()
    plt.legend()
    plt.savefig(extrapolation_path + "/"+match+".png", bbox_inches='tight', dpi=300)

def spectrum_with_extrapolation_execute():
    start_time = time.time()
    shutil.rmtree(spectrum_with_extrapolation_path, ignore_errors=True)
    os.makedirs(spectrum_with_extrapolation_path,exist_ok=True) 
    print_log("--- Spectrum with extrapolation: executing simulations", path=spectrum_with_extrapolation_path)
    mdata = np.genfromtxt(mp_file,delimiter=',')
    print_log("From " + str(mdata[0][0]) + " to " + str(mdata[-1][0]) + " eV", path=spectrum_with_extrapolation_path)
    cmdlineargs = []
    for i in mdata:
        ev = float(i[0])
        lam = 1239.8419842361123824 / (ev * mh)
        mre = float(i[1])/mh
        mim = float(i[2])/mh
        m_abs = math.sqrt(mre**2 + mim**2)
        y_min = (2*math.pi/lam)*(size/grid)*m_abs #y = k*d*|m|
        y_max = 4*y_min
        ys = np.exp(np.linspace(math.log(y_min), math.log(y_max), 9))
        grids = np.rint((2*math.pi/lam)*(size/ys)*m_abs).astype(np.int32)


        os.mkdir(spectrum_with_extrapolation_path + "/" + str(ev))
        with open(spectrum_with_extrapolation_path + "/" + str(ev) + "/ev.csv", 'w') as file:
            csv.writer(file, delimiter=',').writerow([ev])
        for grid_i in grids:
            runpath_i = spectrum_with_extrapolation_path + "/" + str(ev) + "/" + str(grid_i)
            cmdlineargs.append((ev,mre,mim,beam,runpath_i))
    multiprocessing.Pool(procsnumber).starmap(exec_cmdline, cmdlineargs, 1)
    print_log("--- %s seconds" % round((time.time() - start_time),2), path=spectrum_with_extrapolation_path)

def spectrum_with_extrapolation_collect(match):
    print("--- Spectrum with extrapolation: collecting results for", match)
    start_time = time.time()
    global extrapolation_path #this function uses existing extrapolation routine
    temp = extrapolation_path #saving to set this variable back to its value at the end
    evs = sorted([float(d.name) for d in os.scandir(spectrum_with_extrapolation_path) if d.is_dir()])
    values = []
    errors = []
    for ev in evs:
        extrapolation_path = spectrum_with_extrapolation_path + "/" + str(ev)
        extrapolation_collect(match,silent=True)
        data = np.genfromtxt(extrapolation_path + "/"+match+"_fit.csv",delimiter=',')[1:]
        values.append(data[0][0])
        errors.append(data[0][1])
    with open(spectrum_with_extrapolation_path + "/"+match+"_fit.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(["ev",match,"error"])
        writer.writerows(zip(evs,values,errors))
        print("Saved to " + spectrum_with_extrapolation_path + "/"+match+"_fit.csv")
    extrapolation_path = temp #setting variable back to its value
    print("--- %s seconds" % round((time.time() - start_time),2))

def spectrum_with_extrapolation_plot(match):
    #Add exact Mie solution
    # miedata = np.genfromtxt("Peels_mie1545.csv",delimiter=',')
    # plt.plot(miedata[:,0], miedata[:,1], label="Peels_Mie", marker="o", markersize=3, color="red")
    
    print("--- Spectrum with extrapolation: plotting results for", match)
    data = np.genfromtxt(spectrum_with_extrapolation_path + "/"+match+"_fit.csv",delimiter=',')[1:]
    plt.figure()
    plt.plot(data[:,0],data[:,1], label=label_for_plot(match), color="black", marker="x", markersize=5)
    plt.fill_between(data[:,0], data[:,1]-data[:,2], data[:,1]+data[:,2], label="error bar", color="blue", alpha=0.2)
    plt.xlabel("eV")
    plt.show()
    plt.legend()    
    plt.savefig(spectrum_with_extrapolation_path + "/"+match+"_fit.png", bbox_inches='tight', dpi=300)

def scan_execute(ev, step=1):
    start_time = time.time()
    shutil.rmtree(scan_path, ignore_errors=True)
    os.makedirs(scan_path, exist_ok=True)
    print_log("--- Scan: executing simulations", path=scan_path)
    mdata = np.genfromtxt(mp_file,delimiter=',')
    flag = 0
    for line in mdata:
        if float(line[0]) == float(ev):
            flag = 1
            break
    if flag == 0:
        print("ERROR: ev =",ev,"not found in",mp_file)
        return
    mre = float(line[1])/mh
    mim = float(line[2])/mh
    print_log("ev = " + str(ev), path=scan_path)
    print_log("mp_re = " + str(line[1]), path=scan_path)
    print_log("mp_im = " + str(line[2]), path=scan_path)
    
    #adjusting area so the points are exactly in the middle between the dipoles
    d = size/grid #nm
    print_log("dipole size = " + str(d) + " nm", path=scan_path)
    odd = 0.5*(grid % 2)
    left = math.floor(x_left/d + odd) - odd
    right = math.ceil(x_right/d + odd) - odd
    print_log("(x_left,x_right) = (" + str(left) + "," + str(right) + ") dipole sizes", path=scan_path)
    bottom = math.floor(y_bottom/d + odd) - odd
    top = math.ceil(y_top/d + odd) - odd
    print_log("(y_bottom,y_top) = (" + str(bottom) + "," + str(top) + ") dipole sizes", path=scan_path)
    #Step means go through area with 'step' dipole sizes per each step and skip rest of the points.
    x0s = np.linspace(left*d, right*d, round(right - left + 1))[0::step]
    y0s = np.linspace(bottom*d, top*d, round(top - bottom + 1))[0::step]
    # print(x0s)
    # print(y0s)
    
    cmdlineargs = []
    for x0_i in x0s:
        for y0_i in y0s:
            runpath_i = scan_path + "/" + str(x0_i) + "_" + str(y0_i)
            beam = "electron " + str(e_energy) + " " + str(x0_i) + " " + str(y0_i) + " 0 " + str(mh)
            cmdlineargs.append((ev,mre,mim,beam,runpath_i))
    multiprocessing.Pool(procsnumber).starmap(exec_cmdline, cmdlineargs, 1)
    print_log("--- %s seconds" % round((time.time() - start_time),2), path=scan_path)

def scan_collect(match):
    print("--- Scan: collecting results for", match)
    start_time = time.time()
    dirs = sorted([d.name for d in os.scandir(scan_path) if d.is_dir()])
    xs = []
    ys = []
    values = []
    for dir in dirs:
        xy = dir.split("_")
        xs.append(float(xy[0]))
        ys.append(float(xy[1]))
        with open(scan_path + "/" + dir + "/CrossSec-Y", "r") as file:
            for line in file:
                if match in line:
                    value = float(re.findall("[ \t][-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?",line)[0])
                    values.append(value)
                    break
    with open(scan_path + "/"+match+".csv", 'w') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(["x","y",match])
            writer.writerows(zip(xs,ys,values))
            print("Saved to " + scan_path + "/"+match+".csv")
    print("--- %s seconds" % round((time.time() - start_time),2))

def scan_plot(match):
    print("--- Scan: plotting results for", match)
    data = np.genfromtxt(scan_path + "/"+match+".csv",delimiter=',')[1:]
    xs = data[:,0]
    ys = data[:,1]
    zs = data[:,2]
    Nx = len(np.unique(xs))
    Ny = len(np.unique(ys))
    ind = np.lexsort((xs,ys))
    (x,y,z) = (xs[ind].reshape((Ny, Nx)),ys[ind].reshape((Ny, Nx)),zs[ind].reshape((Ny, Nx)))
    # print(x)
    # print(y)
    # print(z)
    plt.figure()
    plt.axes().set_aspect('equal')
    d = size/grid
    plt.imshow(z, extent=(min(xs)-d/2, max(xs)+d/2, min(ys)-d/2, max(ys)+d/2), origin="lower")
    # plt.scatter(x, y, c=z, marker="s") # scatter is the most stable function for visualization, so use this for debugging
    plt.xlabel("x$_0$, nm")
    plt.ylabel("y$_0$, nm")
    plt.colorbar().set_label(label_for_plot(match))
    plt.show()
    plt.savefig(scan_path + "/"+match+".png", bbox_inches='tight', dpi=300)

def label_for_plot(match):
    if match[0] == 'P':
        return match + ", eV$^{-1}$"
    elif match[0] == 'C':
        return match + ", nm$^2$"
    elif match[0] == 'Q':
        return match
    else:
        return match

def print_log(string, path=False):
    print(string)
    if path != False:
        with open(path + "/log.txt", 'a') as file:
            file.write(string + "\n")

# beam = "plane"
# spectrum_execute()
# spectrum_collect("Cext")
# spectrum_plot("Cext")
            
# spectrum_execute()
# spectrum_collect("Peels")
# spectrum_plot("Peels")

# extrapolation_execute(ev=1.9)
# extrapolation_collect("Peels")
# extrapolation_plot("Peels")

spectrum_with_extrapolation_path = "spectrum_with_extrapolation_grid64"
grid = 64
spectrum_with_extrapolation_execute()
spectrum_with_extrapolation_collect("Peels")
spectrum_with_extrapolation_plot("Peels")

# scan_execute(ev=1.45, step=1)
# scan_collect("Peels")
# scan_plot("Peels")



















