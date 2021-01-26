import os, shutil, re, csv, time, multiprocessing, tqdm, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def label_for_plot(match):
    if match[0] == 'P':
        return match + ", eV$^{-1}$"
    elif match[0] == 'C':
        return match + ", nm$^2$"
    elif match[0] == 'Q':
        return match
    else:
        return match

def print_log(string, dirname=False):
    print(string)
    if dirname != False:
        with open(dirname + "/log.txt", 'a') as file:
            file.write(string + "\n")

def mp_range_read(mp_file,ev_min,ev_max):
    mdata = np.genfromtxt(mp_file,delimiter=',')
    mdata_slice = mdata[(mdata[:,0] <= ev_min),0]
    ev_min_nearest = max(mdata_slice) if len(mdata_slice) != 0 else mdata[0,0]
    mdata_slice = mdata[(mdata[:,0] >= ev_max),0]
    ev_max_nearest = min(mdata_slice) if len(mdata_slice) != 0 else mdata[-1,0]
    mdata = mdata[(mdata[:,0] >= ev_min_nearest) & (mdata[:,0] <= ev_max_nearest),:]
    return mdata

def mp_single_read(mp_file,ev):
    mdata = np.genfromtxt(mp_file,delimiter=',')
    mdata = min(mdata, key=lambda x: abs(x[0] - ev))
    return mdata

def cmdline_construct(aw_parameters,adda_cmdlineargs):
    cmdline = aw_parameters["adda_exec"]
    for arg in adda_cmdlineargs:
        if arg != "mh":
            cmdline += f" -{arg} {adda_cmdlineargs[arg]}"
    return cmdline

def ev_to_nm(ev,mh):
    return 1239.8419842361123824 / (ev * mh)

def parse_value(file,match):
    with open(file, "r") as file:
            for line in file:
                if match in line:
                    value = float(re.findall("[ \t][-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?",line)[0])
                    break
    return value

def exec_cmdlines(cmdlines,parallel_procs):
    pool = multiprocessing.Pool(parallel_procs)
    result_list_tqdm = []
    for result in tqdm.tqdm(pool.imap(os.system, cmdlines, 1), total=len(cmdlines)):
        result_list_tqdm.append(result)
    return result_list_tqdm

def spectrum_execute(aw_parameters,adda_cmdlineargs,dirname):
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    print_log("--- Spectrum: executing simulations",dirname)
    if "lambda" in adda_cmdlineargs:
        del adda_cmdlineargs["lambda"]
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    mp_file = aw_parameters["mp_file"]
    ev_min, ev_max = aw_parameters["ev_range"]
    mh = adda_cmdlineargs["mh"]
    mdata = mp_range_read(mp_file,ev_min,ev_max)
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_file: {mp_file}",dirname)
    print_log(f"Varying energy from {mdata[0][0]} to {mdata[-1][0]} eV",dirname)
    cmdlines = []
    for i in mdata:
        cmdline_i = cmdline
        cmdline_i += f" -dir {dirname}/{i[0]}"
        cmdline_i += " -lambda %s" % ev_to_nm(i[0],mh)
        cmdline_i += f" -m {i[1]/mh} {i[2]/mh}"
        cmdline_i += " > /dev/null"
        cmdlines.append(cmdline_i)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)
    
def spectrum_collect(match,dirname):
    evs = sorted([float(d.name) for d in os.scandir(dirname) if d.is_dir()])
    values = []
    for ev in evs:
        values.append(parse_value(f"{dirname}/{ev}/CrossSec-Y",match))
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(["ev",match])
        writer.writerows(zip(evs,values))
        print_log(f"Saved {dirname}/{match}.csv")

def spectrum_plot(match,dirname):
    data = np.genfromtxt(f"{dirname}/{match}.csv", delimiter=',')[1:]
    plt.ion()
    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(data[:,0], data[:,1], label=label_for_plot(match))
    ax.set_xlabel("eV")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(which="both", axis="x", linestyle="dotted")
    ax.legend()
    plt.draw()
    plt.pause(0.001)
    fig.savefig(f"{dirname}/{match}.pdf", bbox_inches='tight')
    print_log(f"Saved {dirname}/{match}.pdf")

def extrapolation_execute(aw_parameters,adda_cmdlineargs,dirname):
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    print_log("--- Extrapolation: executing simulations",dirname)
    mp_file = aw_parameters["mp_file"]
    ev = aw_parameters["ev"]
    mdata = mp_single_read(mp_file,ev)
    mh = adda_cmdlineargs["mh"]
    lam = ev_to_nm(ev,mh)
    adda_cmdlineargs["lambda"] = lam
    mre = mdata[1]/mh
    mim = mdata[2]/mh
    adda_cmdlineargs["m"] = f"{mre} {mim}"
    m_abs = math.sqrt(mre**2 + mim**2)
    size = adda_cmdlineargs["size"]
    grid = adda_cmdlineargs["grid"]
    y_min = (2*math.pi/lam)*(size/grid)*m_abs #y = k*d*|m|
    y_max = 4*y_min
    ys = np.exp(np.linspace(math.log(y_min), math.log(y_max), 9))
    grids = np.rint((2*math.pi/lam)*(size/ys)*m_abs).astype(np.int32)
    del adda_cmdlineargs["grid"]
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_file: {mp_file}",dirname)
    print_log(f"ev = {mdata[0]}",dirname)
    print_log(f"mp_re = {mdata[1]}",dirname)
    print_log(f"mp_im = {mdata[2]}",dirname)
    print_log(f"Varying grids: {grids}",dirname)
    cmdlines = []
    for i in grids:
        cmdline_i = cmdline
        cmdline_i += f" -dir {dirname}/{i}"
        cmdline_i += f" -grid {i}"
        cmdline_i += " > /dev/null"
        cmdlines.append(cmdline_i)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    with open(f"{dirname}/adda_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file)
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)

def extrapolation_collect(match, dirname):
    with open(f"{dirname}/adda_cmdlineargs.csv") as file:
        reader = csv.reader(file)
        adda_cmdlineargs = dict(reader)
    mre = float(adda_cmdlineargs["m"].split(" ")[0])
    mim = float(adda_cmdlineargs["m"].split(" ")[1])
    m_abs = math.sqrt(mre**2 + mim**2)
    grids = np.array(sorted([int(d.name) for d in os.scandir(dirname) if d.is_dir()]))
    values = []
    for grid_i in grids:
        values.append(parse_value(f"{dirname}/{grid_i}/CrossSec-Y",match))
    ys = (2*math.pi/float(adda_cmdlineargs["lambda"]))*(float(adda_cmdlineargs["size"])/grids)*m_abs #y = k*d*|m|
    weights = ys**-3
    fit,cov = np.polyfit(ys, values, 2, w=weights, cov=True)
    a = np.flip(fit)
    error = 2*np.sqrt(np.flip(np.diag(cov)))
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(["grids","ys","values"])
        writer.writerows(zip(grids,ys,values))
        print(f"Saved to {dirname}/{match}.csv")
    with open(f"{dirname}/{match}_fit.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(["a[i]","error[i]"])
        writer.writerows(zip(a,error))
        print(f"Saved to {dirname}/{match}_fit.csv")
    
def extrapolation_plot(match, dirname):
    data = np.genfromtxt(f"{dirname}/{match}.csv",delimiter=',')[1:]
    plt.ion()
    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(data[:,1], data[:,2], label=label_for_plot(match), marker="o", linestyle="none")
    ys_fitted = np.linspace(data[:,1][0],0,100)
    results_fit = np.genfromtxt(f"{dirname}/{match}_fit.csv",delimiter=',')[1:]
    a = results_fit[:,0]
    error = results_fit[:,1]
    points_fitted = a[0] + a[1]*ys_fitted + a[2]*ys_fitted**2
    ax.plot(ys_fitted, points_fitted, label=label_for_plot(match), color="black")
    ax.errorbar(0, a[0], yerr=error[0], color="black", linestyle="", marker="s", capsize=3, barsabove=True, label = "Error bar")
    ax.set_xlabel("y = kd|m|")
    ax.legend()
    plt.draw()
    plt.pause(0.001)
    plt.savefig(f"{dirname}/{match}.pdf", bbox_inches='tight')
    print_log(f"Saved {dirname}/{match}.pdf")

def spectrum_with_extrapolation_execute(aw_parameters,adda_cmdlineargs,dirname):
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname,exist_ok=True) 
    print_log("--- Spectrum with extrapolation: executing simulations",dirname)
    mp_file = aw_parameters["mp_file"]
    ev_min, ev_max = aw_parameters["ev_range"]
    mdata = mp_range_read(mp_file,ev_min,ev_max)
    mh = adda_cmdlineargs["mh"]
    size = adda_cmdlineargs["size"]
    grid = adda_cmdlineargs["grid"]
    ev = float(mdata[0][0])
    lam = ev_to_nm(ev,mh)
    mre = float(mdata[0][1])/mh
    mim = float(mdata[0][2])/mh
    m_abs = math.sqrt(mre**2 + mim**2)
    y_min = (2*math.pi/lam)*(size/grid)*m_abs #y = k*d*|m|
    y_max = 4*y_min
    ys = np.exp(np.linspace(math.log(y_min), math.log(y_max), 9))
    grids = np.rint((2*math.pi/lam)*(size/ys)*m_abs).astype(np.int32) #using same grids for all ev - they would not change
    adda_cmdlineargs["lambda"] = lam
    adda_cmdlineargs["m"] = f"{mre} {mim}"
    with open(f"{dirname}/adda_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file)
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    del adda_cmdlineargs["lambda"]
    del adda_cmdlineargs["m"]
    del adda_cmdlineargs["grid"]
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_file: {mp_file}",dirname)
    print_log(f"Varying energy from {mdata[0][0]} to {mdata[-1][0]} eV",dirname)
    print_log(f"Varying grids: {grids}",dirname)
    cmdlines = []
    for grid_i in grids:
        os.mkdir(f"{dirname}/{grid_i}")
        for mdata_j in mdata:
            cmdline_i = cmdline
            cmdline_i += f" -dir {dirname}/{grid_i}/{mdata_j[0]}"
            cmdline_i += f" -grid {grid_i}"
            cmdline_i += " -lambda %s" % ev_to_nm(mdata_j[0],mh)
            cmdline_i += f" -m {mdata_j[1]/mh} {mdata_j[2]/mh}"
            cmdline_i += " > /dev/null"
            cmdlines.append(cmdline_i)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)

def spectrum_with_extrapolation_collect(match, dirname):
    with open(f"{dirname}/adda_cmdlineargs.csv") as file:
        reader = csv.reader(file)
        adda_cmdlineargs = dict(reader)
    mre = float(adda_cmdlineargs["m"].split(" ")[0])
    mim = float(adda_cmdlineargs["m"].split(" ")[1])
    m_abs = math.sqrt(mre**2 + mim**2)
    grids = np.array(sorted([int(d.name) for d in os.scandir(dirname) if d.is_dir()]))
    evs = sorted([float(d.name) for d in os.scandir(f"{dirname}/{grids[0]}") if d.is_dir()])
    ys = (2*math.pi/float(adda_cmdlineargs["lambda"]))*(float(adda_cmdlineargs["size"])/grids)*m_abs #y = k*d*|m|
    weights = ys**-3
    fit_values = []
    fit_errors = []
    for ev_i in evs:
        values = []
        for grid_j in grids:
            values.append(parse_value(f"{dirname}/{grid_j}/{ev_i}/CrossSec-Y",match))
        fit,cov = np.polyfit(ys, values, 2, w=weights, cov=True)
        a = np.flip(fit)
        error = 2*np.sqrt(np.flip(np.diag(cov)))
        fit_values.append(a[0])
        fit_errors.append(error[0])
    with open(f"{dirname}/{match}_fit.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(["ev",match,"error"])
        writer.writerows(zip(evs,fit_values,fit_errors))
        print(f"Saved to {dirname}/{match}_fit.csv")

def spectrum_with_extrapolation_plot(match,dirname):
    # #Add exact Mie solution
    # miedata = np.genfromtxt(f"Peels_mie.csv",delimiter=',')
    # plt.plot(miedata[:,0], miedata[:,1], label="Peels_Mie", marker="o", markersize=3, color="red")

    data = np.genfromtxt(f"{dirname}/{match}_fit.csv",delimiter=',')[1:]
    plt.ion()
    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(data[:,0], data[:,1], label=label_for_plot(match), color="black")
    ax.fill_between(data[:,0], data[:,1]-data[:,2], data[:,1]+data[:,2], label="error bar", color="blue", alpha=0.2)
    ax.set_xlabel("eV")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(which="both", axis="x", linestyle="dotted")
    ax.legend()  
    plt.draw()
    #plt.pause(0.001)
    plt.savefig(f"{dirname}/{match}_fit.pdf", bbox_inches='tight')
    print_log(f"Saved {dirname}/{match}_fit.pdf")

def scan_execute(aw_parameters,adda_cmdlineargs,dirname):
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    print_log("--- Scan: executing simulations",dirname)
    mp_file = aw_parameters["mp_file"]
    ev = aw_parameters["ev"]
    mdata = mp_single_read(mp_file,ev)
    mh = adda_cmdlineargs["mh"]
    mre = mdata[1]/mh
    mim = mdata[2]/mh
    adda_cmdlineargs["lambda"] = ev_to_nm(ev,mh)
    adda_cmdlineargs["m"] = f"{mre} {mim}"
    size = adda_cmdlineargs["size"]
    grid = adda_cmdlineargs["grid"]
    x_left,x_right = aw_parameters["scan_x_range"]
    y_bottom,y_top = aw_parameters["scan_y_range"]
    step = aw_parameters["scan_step"]
    #adjusting area so the points are exactly in the middle between the dipoles
    d = size/grid #nm
    print_log(f"dipole size = {d} nm",dirname)
    odd = 0.5*(grid % 2)
    left = math.floor(x_left/d + odd) - odd
    right = math.ceil(x_right/d + odd) - odd
    bottom = math.floor(y_bottom/d + odd) - odd
    top = math.ceil(y_top/d + odd) - odd
    x0s = np.linspace(left*d, right*d, round(right - left + 1))[0::step]
    y0s = np.linspace(bottom*d, top*d, round(top - bottom + 1))[0::step]
    # print(x0s)
    # print(y0s)
    beam_list = adda_cmdlineargs["beam"].split(" ")
    del adda_cmdlineargs["beam"]
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    with open(f"{dirname}/adda_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file)
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_file: {mp_file}",dirname)
    print_log(f"ev = {ev}",dirname)
    print_log(f"mp_re = {mdata[1]}",dirname)
    print_log(f"mp_im = {mdata[2]}",dirname)
    print_log(f"Varying (x_left,x_right) = ({left},{right}) dipole sizes",dirname)
    print_log(f"Varying (y_bottom,y_top) = ({bottom},{top}) dipole sizes",dirname)
    cmdlines = []
    for x0_i in x0s:
        for y0_i in y0s:
            cmdline_i = cmdline
            cmdline_i += f" -dir {dirname}/{x0_i}_{y0_i}"
            beam_list[2], beam_list[3] = str(x0_i), str(y0_i)
            beam = (" ").join(beam_list)
            cmdline_i += f" -beam {beam}"
            cmdline_i += " > /dev/null"
            cmdlines.append(cmdline_i)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)

def scan_collect(match, dirname):
    dirs = sorted([d.name for d in os.scandir(dirname) if d.is_dir()])
    xs = []
    ys = []
    values = []
    for dir_i in dirs:
        xy = dir_i.split("_")
        xs.append(float(xy[0]))
        ys.append(float(xy[1]))
        values.append(parse_value(f"{dirname}/{dir_i}/CrossSec-Y",match))
    with open(f"{dirname}/{match}.csv", 'w') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(["x","y",match])
            writer.writerows(zip(xs,ys,values))
            print(f"Saved to {dirname}/{match}.csv")

def scan_plot(match,dirname):
    with open(f"{dirname}/adda_cmdlineargs.csv") as file:
        reader = csv.reader(file)
        adda_cmdlineargs = dict(reader)
    size = float(adda_cmdlineargs["size"])
    grid = float(adda_cmdlineargs["grid"])
    data = np.genfromtxt(f"{dirname}/{match}.csv",delimiter=',')[1:]
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
    plt.ion()
    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect('equal')
    d = size/grid
    plt.imshow(z, extent=(min(xs)-d/2, max(xs)+d/2, min(ys)-d/2, max(ys)+d/2), origin="lower", cmap="rainbow")
    # plt.scatter(x, y, c=z, marker="s") # scatter is the most stable function for visualization, so use this for debugging
    ax.set_xlabel("x$_0$, nm")
    ax.set_ylabel("y$_0$, nm")
    plt.colorbar().set_label(label_for_plot(match))
    #plt.axis('off')
    plt.draw()
    #plt.pause(0.001)
    plt.savefig(f"{dirname}/{match}.pdf", bbox_inches='tight')
    print_log(f"Saved {dirname}/{match}_fit.pdf")

