# run v0.2
import os, shutil, re, csv, time, multiprocessing, math
import numpy as np
import matplotlib.pyplot as plt

addaexec = "/Users/user/Documents/GitHub/adda/src/seq/adda" #path to adda executable
runpath = "/Users/user/Documents/e_field/spectrum" #where runs are stored
q = -4.803204673e-10
c = 29979245800
rest = 510.99895
hbar = 1.054571817e-27
hbar_ev = 6.582119569e-16

# Particle
shape = "sphere"
#shape = "prism 3 0.1282"
size = 10 #nm
grid = 16
mfile = "/Users/user/Documents/e_field/spectrum/Ag_Palik_Garcia330-55.csv" #m_particle: each string contains: ev n k
mh = 1.5 #m_host
mh1 = mh

# EELS and CL
r0 = (6,0,0) #nm
#prop = (0,1,1) #propagation direction vector
elenergy = 100 #keV
beam = "electron " + str(elenergy) + " " + str(r0[0]) + " " + str(r0[1]) + " " + str(r0[2]) + " " + str(mh)

# Precision and performance
eps = 5 #Residual
calc_extrapolation = False
procsnumber = multiprocessing.cpu_count() #number of parallel processes is equal to the number of processor cores
#procsnumber = 2 #manually select number of processes

def exec_cmdline(ev,mpre,mpim,beam,runpath_i,matches):
    lam = 1239.8419842361123824 / (ev * mh)
    #lam = 1239.8419842361123824 / ev
    (mre,mim) = (mpre/mh,mpim/mh)
    cmdline = addaexec + \
    " -shape " + shape + \
    " -size " + str(size) + \
    " -lambda " + str(lam) + \
    " -m " + str(mre) + " " + str(mim) + \
    " -beam " + beam + \
    " -eps " + str(eps) + \
    " -pol cm -sym enf -scat_matr none -no_vol_cor" + \
    " -dir '" + runpath_i + "'"

    cmdline += " -grid " + str(grid)
    cmdline += " > /dev/null"
    flag = os.system(cmdline)
    if flag != 0:
        print("'" + cmdline + "' ran with exit code ", flag)
        return -1
    values = []
    for match in matches:
        match = match + "\t"
        #with open(runpath_i + "/log", "r") as file:
        with open(runpath_i + "/CrossSec-Y", "r") as file:
            for line in file:
                if match in line:
                    value = float(re.findall("[ \t][-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?",line)[0])
                    values.append(value)
                    break
    return values

def writefile(matches,var,result):
    result = np.array(result)
    timestamp = time.strftime('%Y-%m-%d_%H%M%S', time.localtime())
    for index,match in enumerate(matches):
        with open(runpath + "/temp/result_"+match+".csv", 'w') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerows(zip(*var,result[:,index]))
        with open(runpath + "/history/result_"+match+"_"+timestamp+".csv", 'w') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerows(zip(*var,result[:,index]))

def run_spectrum(matches):
    print("--- Running spectrum for",matches)
    start_time = time.time()
    mdata = np.genfromtxt(mfile,delimiter=',')
    evs = mdata[:,0]
    cmdlineargs = []
    for i in mdata:
        runpath_i = runpath + "/output/" + str(i[0])
        cmdlineargs.append((i[0],i[1],i[2],beam,runpath_i,matches))
    shutil.rmtree(runpath + "/output", ignore_errors=True)
    os.mkdir(runpath + "/output") 
    pool = multiprocessing.Pool(procsnumber)
    result = pool.starmap(exec_cmdline, cmdlineargs, 1)
    writefile(matches,(evs,),result)
    print("--- %s seconds" % round((time.time() - start_time),2))

def plot_spectrum(matches):
    for match in matches:
        data = np.genfromtxt(runpath + "/temp/result_"+match+".csv",delimiter=',')
        xs = data[:,0]
        ys = data[:,1]
        plt.plot(xs, ys, label = match)
        plt.show()
    plt.legend()    
    plt.savefig(runpath + "/temp/result_spectrum.png", bbox_inches='tight')

def run_scan(matches):
    print("--- Running scan for",matches)
    start_time = time.time()
    ppnm = 1 #points per nm
    dpnm = 1 #dipoles per nm
    global grid
    grid = size*dpnm
    (x_start,x_stop) = (0,50) #nm
    x_steps = int(abs(x_stop - x_start)*ppnm + 1) #points
    (y_start,y_stop) = (-0,50) #nm
    y_steps = int(abs(y_stop - y_start)*ppnm + 1)
    ev = 3.6
    with open(mfile, "r") as file:
        for line in file:
            if str(ev)+"," in line:
                mpre = float(re.findall("[-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?",line)[1])
                mpim = float(re.findall("[-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?",line)[2])
                break
    x0s = np.linspace(x_start, x_stop, x_steps)
    y0s = np.linspace(y_start, y_stop, y_steps)
    xs = []
    ys = []
    cmdlineargs = []
    for x0_i in x0s:
        for y0_i in y0s:
            runpath_i = runpath + "/output/" + str(x0_i) + "_" + str(y0_i)
            beam = "electron " + str(elenergy) + " " + str(x0_i) + " " + str(y0_i) + " 0 " + str(mh)
            cmdlineargs.append((ev,mpre,mpim,beam,runpath_i,matches))
            shutil.rmtree(runpath_i, ignore_errors=True)
            xs.append(float(x0_i))
            ys.append(float(y0_i))
    shutil.rmtree(runpath + "/output/", ignore_errors=True)
    os.mkdir(runpath + "/output") 
    pool = multiprocessing.Pool(procsnumber)
    result = pool.starmap(exec_cmdline, cmdlineargs)
    writefile(matches,(xs,ys),result)
    print("--- %s seconds" % round((time.time() - start_time),2))

def plot_scan(matches):
    match = matches[0]
    data = np.genfromtxt(runpath + "/temp/result_"+match+".csv",delimiter=',')
    xs = data[:,0]
    ys = data[:,1]
    zs = data[:,2]
    Nx = len(np.unique(xs))
    Ny = len(np.unique(ys))
    x = xs.reshape((Nx, Ny))
    y = ys.reshape((Nx, Ny))
    z = zs.reshape((Nx, Ny))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    #plt.scatter(x, y, c=z)
    plt.pcolormesh(x, y, z)
    plt.colorbar()
    plt.savefig(runpath + "/temp/result_"+match+".png", bbox_inches='tight')

matches = ["Peels"]

#run_spectrum(matches)
#plot_spectrum(matches)

'''
if calc_extrapolation == True:
    m_abs = math.sqrt(mre**2 + mim**2)
    y_min = (2*math.pi/lam)*(size/grid)*m_abs #y = k*d*|m|
    y_max = 4*y_min
    ys = np.exp(np.linspace(math.log(y_min), math.log(y_max), 9))
    grids = np.round((2*math.pi/lam)*(size/ys)*m_abs)
    print(grids)
'''

'''
# lam не делить на mh
size *= mh
r0 = [i * mh for i in r0]
elenergy = 510.99895*(1/sqrt(1-mh*mh*(1-(510.99895/(510.99895 + elenergy))**2)) - 1)
beam = "electron " + str(elenergy) + " " + str(r0[0]) + " " + str(r0[1]) + " " + str(r0[2]) + " 1"
mh1 = 1
print(size, r0, elenergy)
# P отсклейить обратно
run_spectrum(matches)
plot_spectrum(matches)
'''

#run_scan(["Peels"])
#plot_scan(["Peels"])






















