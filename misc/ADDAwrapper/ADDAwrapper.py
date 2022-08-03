import os, shutil, re, csv, time, multiprocessing, tqdm, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter
 
def addaexec_find(mode="seq"):
    if mode!="seq" and mode!="mpi" and mode!="ocl":
        print(f"ERROR: unkwnown mode {mode}")
        return
    if mode == "seq": 
        name="adda"
    else: 
        name=f"adda_{mode}" 
    if os.system(name  + " > " + os.devnull + " 2>&1") == 0 : return name
    cmdline = os.path.abspath(__file__ + f"/../../../src/seq/{name}")
    if os.system(cmdline  + " -V > " + os.devnull + " 2>&1") == 0 : return cmdline
    cmdline = os.path.abspath(__file__ + f"/../../../win64/{name}")
    if os.name == "nt" and os.system(cmdline  + " -V > " + os.devnull + " 2>&1") == 0 : return cmdline
    print(f"ERROR: No working {name} binary found")

def label_for_plot(match):
    if match[0] == "P":
        if match[1:] == "ext": return r"$\it{PCQ}$" + r"$_{\rm ext}$" + ", eV$^{-1}$"
        return r"$\it{" + match[0] + r"}$" + r"$_{\rm " + match[1:].upper() + "}$" + ", eV$^{-1}$"
    elif match[0] == "C":
        return r"$\it{" + match[0] + r"}$" + r"$_{\rm " + match[1:] + "}$" + ", nm$^2$"
    elif match[0] == "Q":
        return r"$\it{" + match[0] + r"}$" + r"$_{\rm " + match[1:] + "}$"
    else:
        return match

def label_for_plot_arbunits(match):
    if match[0] == "P":
        return match[0] + r"$_{\rm " + match[1:].upper() + "}$" + ", arb. units"
    elif match[0] == "C":
        return match + ", arb. units"
    elif match[0] == "Q":
        return match
    else:
        return match

def color_for_plot(match):
    if match == "Peels":
        return "royalblue"
    elif match == "Pcl":
        return "violet"
    else:
        return "black"

def print_log(string, dirname=False, silent=False):
    if silent == False:
        print(string)
    if dirname != False:
        with open(dirname + "/log.txt", 'a') as file:
            file.write(string + "\n")

def delkeys_silent(dictionary, keys):
    for key in keys:
        if key in dictionary:
            del dictionary[key]

def mp_range_read(mp_file,ev_min,ev_max):
    mdata = np.genfromtxt(mp_file,delimiter=',')
    mdata_slice = mdata[(mdata[:,0] <= ev_min),0]
    ev_min_nearest = max(mdata_slice) if len(mdata_slice) != 0 else mdata[0,0]
    mdata_slice = mdata[(mdata[:,0] >= ev_max),0]
    ev_max_nearest = min(mdata_slice) if len(mdata_slice) != 0 else mdata[-1,0]
    mdata = mdata[(mdata[:,0] >= ev_min_nearest) & (mdata[:,0] <= ev_max_nearest),:]
    return mdata

def dipole_middles(point1, point2, d, odd):
    (x_left,y_bottom) = point1
    (x_right,y_top) = point2
    left = math.floor(x_left/d + odd) - odd
    right = math.ceil(x_right/d + odd) - odd
    bottom = math.floor(y_bottom/d + odd) - odd
    top = math.ceil(y_top/d + odd) - odd
    #print(left,right,bottom,top)
    x0s = np.linspace(left*d, right*d, abs(round(right - left))+1)
    x0s = np.around(x0s,8)
    y0s = np.linspace(bottom*d, top*d, abs(round(top - bottom))+1)
    y0s = np.around(y0s,8)
    return x0s, y0s

def mp_single_read(mp_file,ev):
    mdata = np.genfromtxt(mp_file,delimiter=',')
    mdata = min(mdata, key=lambda x: abs(x[0] - ev))
    return mdata

def cmdline_construct(aw_parameters,adda_cmdlineargs):
    cmdline = aw_parameters["adda_exec"]
    for arg in adda_cmdlineargs:
        cmdline += f" -{arg} {adda_cmdlineargs[arg]}"
    return cmdline

def ev_to_nm(ev):
    return 1239.8419842361123824 / ev

def nm_to_ev(nm):
    return 1239.8419842361123824 / nm

def parse_value(file,match):
    with open(file, "r") as file:
            for line in file:
                if match in line:
                    value = float(re.findall("[ \t][-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?",line)[0])
                    break
    return value

def midpoints(x):
    sl = ()
    for i in range(x.ndim):
        x = (x[sl + np.index_exp[:-1]] + x[sl + np.index_exp[1:]]) / 2.0
        sl += np.index_exp[:]
    return x

def geom_gen(aw_parameters,adda_cmdlineargs):
        cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
        #cmdline = aw_parameters["adda_exec"] 
        cmdline += f" -dir {aw_parameters['dirname']}/geom"
        cmdline += " -prognosis"
        #cmdline += " -no_vol_cor"
        #cmdline += f" -shape {adda_cmdlineargs['shape']}"
        #cmdline += f" -grid {adda_cmdlineargs['grid']}"
        cmdline += " -save_geom"
        cmdline += " > " + os.devnull
        os.system(cmdline)
        #print(cmdline)

def geometry(geom_path,colorlist=None,half=False): 
    skip_header = 0
    with open(geom_path, "r") as file:
            for line in file:
                if "Nmat" in line:
                    skip_header=1
    data = np.genfromtxt(geom_path, delimiter=' ', dtype="int", skip_header=skip_header)
    if half == True:
        #data = data[(data[:,0]<.5*(max(data[:,0]) + min(data[:,0]))) | (data[:,2]<.5*(max(data[:,2]) + min(data[:,2])))] #quarter
        data = data[data[:,0]<.5*(max(data[:,0]) + min(data[:,0]))] #half
    xs, ys, zs = data[:,0], data[:,1], data[:,2]
    xs, ys, zs = xs-min(xs), ys-min(ys), zs-min(zs)
    x, y, z = np.indices((max(xs)-min(xs)+2,max(ys)-min(ys)+2,max(zs)-min(zs)+2)).astype(float)
    xc = midpoints(x)
    voxels = np.zeros(xc.shape, dtype=bool)
    colors = np.empty(xc.shape, dtype=object)
    if data.shape[1] == 4:
        ms = data[:,3] - 1
    else:
        if colorlist == None:
            ms = np.ones(xs.shape, dtype="int")
        else:
            ms = np.zeros(xs.shape, dtype="int")
    if colorlist == None:
        colorlist = ["deepskyblue", "silver", "gold", "yellowgreen", "tomato", "darkviolet", "peru", "lime"]
    for i in range(len(xs)):
        voxels[xs[i],ys[i],zs[i]] = True
        colors[xs[i],ys[i],zs[i]] = colorlist[ms[i]]
            
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    #print(voxels)
    ax.voxels(x, y, z, voxels, facecolors=colors, edgecolors="grey", linewidth=0.05, alpha=0.9, zorder=1)
    
    ax.set_xlim([0, xc.shape[0]])
    ax.set_ylim([0, xc.shape[1]])
    ax.set_zlim([0, xc.shape[2]])
    ax.set(xlabel='x', ylabel='y', zlabel='z')
    ax.set_box_aspect([max(xs)-min(xs)+1,max(ys)-min(ys)+1,max(zs)-min(zs)+1])
    ax.view_init(30, -120)
    return fig,ax

def plot_create():
    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    #ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0,0))
    ax.yaxis.set_major_formatter(formatter)
    ax.set_xlabel("Energy, eV")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(which="both", axis="x", linestyle="dotted", zorder=0)
    ax.tick_params(bottom=True, top=True, left=True, right=True, which = "both")
    ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    plot_setrcparams()
    return fig, ax

def plot_setaspect(ax):
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ytop-ybottom))*.5625)

def plot_setrcparams():
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 16
    
    # SMALL_SIZE = 16
    # MEDIUM_SIZE = 18
    # BIGGER_SIZE = 16
    
    # SMALL_SIZE = 24
    # MEDIUM_SIZE = 32
    # BIGGER_SIZE = 32
    
    plt.rc('font', **{'family': 'serif', 'serif': 'Arial'})
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    
    #plt.rc('legend', fontsize=32)    # legend fontsize

def exec_cmdlines(cmdlines,parallel_procs):
    pool = multiprocessing.Pool(parallel_procs)
    result_list_tqdm = []
    for result in tqdm.tqdm(pool.imap(os.system, cmdlines, 1), total=len(cmdlines)):
        result_list_tqdm.append(result)
    return result_list_tqdm

def varyany_execute(aw_parameters,adda_cmdlineargs,dirname,var,var_range):
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    os.makedirs(f"{dirname}/ADDA_output", exist_ok=True)
    print()
    print_log("--- Vary any parameter: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    if var in adda_cmdlineargs:
        del adda_cmdlineargs[var]
    mp_file = aw_parameters["mp_file"]
    ev = aw_parameters["ev"]
    mdata = mp_single_read(mp_file,ev)
    lam = ev_to_nm(ev)
    adda_cmdlineargs["lambda"] = lam
    mre = mdata[1]
    mim = mdata[2]
    adda_cmdlineargs["m"] = f"{mre} {mim}"
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    print_log(f"{cmdline}",dirname)
    print_log(f"mp = {mre} + {mim}*I",dirname)
    print_log(f"Varying {var} from {var_range[0]} to {var_range[-1]}",dirname)
    cmdlines = []
    for i in var_range:
        cmdline_i = cmdline
        output_dir = os.path.abspath(dirname + f"/ADDA_output/{i}")
        cmdline_i += f' -dir "{output_dir}"'
        cmdline_i += f" -{var} %s" % str(i)
        cmdline_i += " > " + os.devnull
        cmdlines.append(cmdline_i)
    #print(cmdlines)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)
    
def varyany_collect(match,dirname, silent=False):
    xs = sorted([d.name for d in os.scandir(f"{dirname}/ADDA_output") if d.is_dir()])
    #print(xs)
    values = []
    for x in xs:
        values.append(parse_value(f"{dirname}/ADDA_output/{x}/CrossSec-Y",match))
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["var",match])
        writer.writerows(zip(xs,values))
        print_log(f"Saved {dirname}/{match}.csv", silent=silent)
        
def varyany_plot(match,dirname):
    data = np.genfromtxt(f"{dirname}/{match}.csv", delimiter=',')[1:]
    fig,ax = plot_create()
    ax.plot(data[:,0], data[:,1], label=label_for_plot(match), color=color_for_plot(match), linestyle="none", marker=".")
    ax.set_xlim([min(data[:,0]),max(data[:,0])])
    ax.set_ylabel(label_for_plot(match))
    #ax.legend()
    plot_setaspect(ax)
    ax.set_xlabel("")
    fig.savefig(f"{dirname}/{match}.svg", bbox_inches='tight')
    print_log(f"Saved {dirname}/{match}.svg")

def spectrum_execute(aw_parameters,adda_cmdlineargs):
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    dirname = aw_parameters["dirname"]
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    os.makedirs(f"{dirname}/ADDA_output", exist_ok=True)
    print()
    print_log("--- Spectrum: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    delkeys_silent(adda_cmdlineargs, ["lambda"])
    if 'ms_files' in aw_parameters.keys():
        surfheight = adda_cmdlineargs["surf"].split(" ")[0]
        delkeys_silent(adda_cmdlineargs, ["surf"])
        msdata = np.genfromtxt(aw_parameters["ms_files"][0],delimiter=',')
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    mp_files = aw_parameters["mp_files"]
    if 'ev_range' in aw_parameters.keys():
        evs = [round(x,2) for x in aw_parameters["ev_range"]]
    else:
        evs = [round(x,2) for x in aw_parameters["nm_range"]]
    mdata = dict()
    for i in range(len(mp_files)):
        mdata[i] = np.genfromtxt(mp_files[i],delimiter=',')
    mps = dict()
    for ev in evs:
        mps[ev] = ""
        for i in mdata:
            nk = min(mdata[i], key=lambda x: abs(x[0] - ev))
            nk = str(nk[1]) + " " + str(nk[2]) + " "
            mps[ev] += str(nk)
    with open(f"{dirname}/ADDA_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_files: {mp_files}",dirname)
    if 'ev_range' in aw_parameters.keys():
        print_log(f"Varying energy from {evs[0]} to {evs[-1]} eV",dirname)
    else:
        print_log(f"Varying wavelength from {evs[0]} to {evs[-1]} nm",dirname)
    cmdlines = []
    for ev in evs:
        cmdline_i = cmdline
        output_dir = os.path.abspath(dirname + f"/ADDA_output/{float(ev)}")
        cmdline_i += f' -dir "{output_dir}"'
        if 'ev_range' in aw_parameters.keys():
            cmdline_i += " -lambda %s" % ev_to_nm(ev)
        else:
            cmdline_i += " -lambda %s" % ev
        cmdline_i += " -m " + f" {mps[ev]} "
        if 'ms_files' in aw_parameters.keys():
            nk = min(msdata, key=lambda x: abs(x[0] - ev))
            nk = str(nk[1]) + " " + str(nk[2])
            cmdline_i += " -surf " + f" {surfheight} {nk}"
        cmdline_i += " > " + os.devnull
        cmdlines.append(cmdline_i)
    #print(cmdlines)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)
    
def spectrum_collect(match,aw_parameters,silent=False):
    dirname = aw_parameters["dirname"]
    evs = sorted([float(d.name) for d in os.scandir(f"{dirname}/ADDA_output") if d.is_dir()])
    values = []
    for ev in evs:
        values.append(parse_value(f"{dirname}/ADDA_output/{ev}/CrossSec-Y",match))
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        if 'ev_range' in aw_parameters.keys():
            writer.writerow(["eV",match])
        else:
            writer.writerow(["nm",match])
        writer.writerows(zip(evs,values))
        print_log(f"Saved {dirname}/{match}.csv", silent=silent)

def spectrum_plot(match,aw_parameters):
    dirname = aw_parameters["dirname"]
    data = np.genfromtxt(f"{dirname}/{match}.csv", delimiter=',')[1:]
    fig,ax = plot_create()
    if 'ev_range' in aw_parameters.keys():
        ax.set_xlabel("Energy, eV")
    else:
        ax.set_xlabel("Wavelength, nm")
    #print(data)
    ax.plot(data[:,0], data[:,1], label=label_for_plot(match), color=color_for_plot(match), linewidth=3)
    ax.set_xlim([min(data[:,0]),max(data[:,0])])
    ax.set_ylabel(label_for_plot(match))
    #ax.legend()
    #plot_setaspect(ax)
    fig.savefig(f"{dirname}/{match}.png", bbox_inches='tight', dpi=600)
    print_log(f"Saved {dirname}/{match}.png")
    return fig,ax
    
def spectrumline_execute(aw_parameters,adda_cmdlineargs):
    dirname = aw_parameters["dirname"]
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    print()
    print_log("--- Spectrum for a set of points on a line: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    mp_files = aw_parameters["mp_files"]
    if 'ev_range' in aw_parameters.keys():
        evs = [round(x,2) for x in aw_parameters["ev_range"]]
    else:
        evs = [round(x,2) for x in aw_parameters["nm_range"]]
    mdata = dict()
    for i in range(len(mp_files)):
        mdata[i] = np.genfromtxt(mp_files[i],delimiter=',')
    
    mps = dict()
    for ev in evs:
        mps[ev] = ""
        for i in mdata:
            nk = min(mdata[i], key=lambda x: abs(x[0] - ev))
            nk = str(nk[1]) + " " + str(nk[2]) + " "
            mps[ev] += str(nk)
    with open(f"{dirname}/ADDA_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    size = adda_cmdlineargs["size"]
    grid = adda_cmdlineargs["grid"]
    point1 = aw_parameters["spectrumline_startpoint"]
    point2 = aw_parameters["spectrumline_endpoint"]
    howmanypoints = aw_parameters["spectrumline_points"]
    x_step, y_step = (point2[0]-point1[0])/(howmanypoints-1), (point2[1]-point1[1])/(howmanypoints-1)
    #adjusting area so the points are exactly in the middle between the dipoles
    d = size/grid #nm
    odd = 0.5*(grid % 2)
    x0s, y0s = dipole_middles(point1,point2,d,odd)
    # print(x0s)
    # print(y0s)
    points = []
    for i in range(0,howmanypoints):
        x0_i = point1[0] + x_step*i
        y0_i = point1[1] + y_step*i
        x_i = min(x0s, key=lambda val: abs(val - x0_i))
        y_i = min(y0s, key=lambda val: abs(val - y0_i))
        #print(i, x0_i, x_i, y0_i, y_i)
        points.append((x_i,y_i))
    points = np.array(points)
    #print(points)
    
    # fig = plt.figure(constrained_layout=True)
    # ax = fig.add_subplot(1, 1, 1)
    # ax.plot(points[:,0], points[:,1], linewidth=3, marker="o")
    # #return
    
    beam_list = adda_cmdlineargs["beam_center"].split(" ")
    delkeys_silent(adda_cmdlineargs, ["lambda","m","beam_center"])
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_files: {mp_files}",dirname)
    print_log(f"dipole size = {d} nm",dirname)
    print_log(f"Varying position from {tuple(points[0,:])} nm to {tuple(points[-1,:])} nm ({howmanypoints} points)",dirname)
    if 'ev_range' in aw_parameters.keys():
        print_log(f"Varying energy from {evs[0]} to {evs[-1]} eV",dirname)
    else:
        print_log(f"Varying wavelength from {evs[0]} to {evs[-1]} nm",dirname)
    
    cmdlines = []
    counter = 0
    for p in points:
        counter += 1
        point_dir = f"{dirname}/" + "{:03d}".format(counter) + f"_{p[0]}_{p[1]}"
        #print(point_dir)
        os.mkdir(point_dir)
        os.mkdir(point_dir+"/ADDA_output")
        beam_list[0], beam_list[1] = str(p[0]), str(p[1])
        beam_center = (" ").join(beam_list)
        #print(beam)
        for ev in evs:
            cmdline_i = cmdline
            cmdline_i += f" -beam_center {beam_center}"
            dir_fixed = os.path.abspath(point_dir + f"/ADDA_output/{float(ev)}")
            cmdline_i += f' -dir "{dir_fixed}"'
            if 'ev_range' in aw_parameters.keys():
                cmdline_i += " -lambda %s" % ev_to_nm(ev)
            else:
                cmdline_i += " -lambda %s" % ev
            cmdline_i += " -m" + f" {mps[ev]} "
            cmdline_i += " > " + os.devnull
            cmdlines.append(cmdline_i)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)
    
def spectrumline_collect(match, aw_parameters):
    dirname = aw_parameters["dirname"]
    dirs = sorted([d.name for d in os.scandir(dirname) if d.is_dir()])
    #print(dirs)
    points = []
    ys = []
    for dir_i in dirs:
        aw_parameters["dirname"] = f"{dirname}/{dir_i}"
        spectrum_collect(match,aw_parameters, silent=True)
        ys_i = np.genfromtxt(f"{dirname}/{dir_i}/{match}.csv", delimiter=',')[1:,1]
        ys.append(ys_i)
        points.append(dir_i.split("_"))
    aw_parameters["dirname"] = dirname
    #points = np.array(points)
    #print(points)
    xs = np.genfromtxt(f"{dirname}/{dirs[0]}/{match}.csv", delimiter=',')[1:,0]
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerows(zip(["Point no.", "x [nm]", "y [nm]"],*points))
        writer.writerow("-"*(len(points)+1))
        if 'ev_range' in aw_parameters.keys():
            valuenames = ["eV"] + [f'{match}']*len(points)
        else:
            valuenames = ["nm"] + [f'{match}']*len(points)
        
        writer.writerow(valuenames)
        writer.writerows(zip(xs,*ys))
        print(f"Saved to {dirname}/{match}.csv")

def spectrumline_plot(match, aw_parameters):
    dirname = aw_parameters["dirname"]
    #alldata = np.genfromtxt(f"{dirname}/{match}.csv", delimiter=',',dtype=None, encoding=None)
    data = np.genfromtxt(f"{dirname}/{match}.csv", delimiter=',')
    fig1,ax1 = plot_create()
    ax1.set_ylabel(label_for_plot_arbunits(match))
    fig2,ax2 = plot_create()
    if 'ev_range' in aw_parameters.keys():
        ax1.set_xlabel("Energy, eV")
        ax2.set_xlabel("Energy, eV")
        ax1.set_xlabel("Энергия, эВ")
        ax2.set_xlabel("Энергия, эВ")
        ax1.set_ylabel(r"P$_{\rm EELS}$" + ", п. е.")
        ax2.set_ylabel(r"P$_{\rm CL}$" + ", п. е.")
    else:
        ax1.set_xlabel("Wavelength, nm")
        ax2.set_xlabel("Wavelength, nm")
        ax1.set_xlabel("Длина волны, нм")
        ax2.set_xlabel("Длина волны, нм")
        ax1.set_ylabel(r"P$_{\rm CL}$" + ", п. е.")
        ax2.set_ylabel(r"P$_{\rm CL}$" + ", п. е.")
    ax2.set_ylabel(label_for_plot(match))
    
    xs = data[5:,0]
    #print(xs)
    ys2 = np.zeros(np.shape(data.T[0][5:]))
    for point in data.T[1:]:
        num = int(point[0])
        ys = point[5:]/max(point[5:]) + .02*num
        ys2 += point[5:]
        #print(num)
        #print(ys)
        ax1.plot(xs, ys, color=color_for_plot(match), linewidth=1.5, zorder=(1-0.001*num))
        ax1.fill_between(xs, min(ys), ys, facecolor="white", alpha=.3, zorder=(1-0.001*num))
    ys2 /= num
    ax1.set_xlim([min(xs),max(xs)])
    ax2.set_xlim([min(xs),max(xs)])
    ax2.plot(xs, ys2, color=color_for_plot(match), linewidth=3)
    #ax1.set_yscale('log')
    plot_setaspect(ax1)
    plot_setaspect(ax2)
    fig1.savefig(f"{dirname}/{match}.png", bbox_inches='tight', dpi=600)
    print_log(f"Saved {dirname}/{match}.png")
    fig2.savefig(f"{dirname}/{match}_averaged.png", bbox_inches='tight', dpi=600)
    print_log(f"Saved {dirname}/{match}_averaged.png")
    return fig1,ax1

def extrapolation_execute(aw_parameters,adda_cmdlineargs):
    dirname = aw_parameters["dirname"]
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    os.makedirs(f"{dirname}/ADDA_output", exist_ok=True)
    print()
    print_log("--- Extrapolation: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    mp_files = aw_parameters["mp_files"]
    ev = aw_parameters["ev"]
    mdata = dict()
    for i in range(len(mp_files)):
        mdata[i] = np.genfromtxt(mp_files[i],delimiter=',')
    mps = dict()
    mps[ev] = ""
    m_abs = 0
    for i in mdata:
        nk = min(mdata[i], key=lambda x: abs(x[0] - ev))
        if m_abs < math.sqrt(float(nk[1])**2 + float(nk[2])**2):
            m_abs = math.sqrt(float(nk[1])**2 + float(nk[2])**2)
        nk = str(nk[1]) + " " + str(nk[2]) + " "
        mps[ev] += str(nk) + " "
        
    lam = ev_to_nm(ev)
    adda_cmdlineargs["lambda"] = lam
    adda_cmdlineargs["m"] = f"{mps[ev]}"
    size = adda_cmdlineargs["size"]
    grid = adda_cmdlineargs["grid"]
    y_min = (2*math.pi/lam)*(size/grid)*m_abs #y = k*d*|m|
    y_max = 4*y_min
    ys = np.exp(np.linspace(math.log(y_min), math.log(y_max), 9))
    grids = np.rint((2*math.pi/lam)*(size/ys)*m_abs).astype(np.int32)
    #grids = np.asarray([128, 112, 96, 80, 64, 56, 48, 40, 32])
    del adda_cmdlineargs["grid"]
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_files: {mp_files}",dirname)
    print_log(f"ev = {ev}",dirname)
    print_log(f"m_p = {mps[ev]}",dirname)
    print_log(f"Varying grids: {grids}",dirname)
    cmdlines = []
    for i in grids:
        cmdline_i = cmdline
        dir_fixed = os.path.abspath(dirname + f"/ADDA_output/{i}")
        cmdline_i += f' -dir "{dir_fixed}"'
        cmdline_i += f" -grid {i}"
        cmdline_i += " > " + os.devnull
        cmdlines.append(cmdline_i)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    with open(f"{dirname}/ADDA_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)

def extrapolation_collect(match, aw_parameters):
    dirname = aw_parameters["dirname"]
    with open(f"{dirname}/ADDA_cmdlineargs.csv") as file:
        reader = csv.reader(file)
        adda_cmdlineargs = dict(reader)
    m_abs = 0
    mdata = re.split(" +",adda_cmdlineargs["m"])
    print(mdata)
    for i in range(int(len(mdata)/2)):
        if m_abs < math.sqrt(float(mdata[2*i])**2 + float(mdata[2*i+1])**2):
            m_abs = math.sqrt(float(mdata[2*i])**2 + float(mdata[2*i+1])**2)
    grids = np.array(sorted([int(d.name) for d in os.scandir(f"{dirname}/ADDA_output") if d.is_dir()]))
    values = []
    for grid_i in grids:
        values.append(parse_value(f"{dirname}/ADDA_output/{grid_i}/CrossSec-Y",match))
    ys = (2*math.pi/float(adda_cmdlineargs["lambda"]))*(float(adda_cmdlineargs["size"])/grids)*m_abs #y = k*d*|m|
    weights = ys**-3
    fit,cov = np.polyfit(ys, values, 2, w=weights, cov=True)
    a = np.flip(fit)
    error = 2*np.sqrt(np.flip(np.diag(cov)))
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["grids","ys","values"])
        writer.writerows(zip(grids,ys,values))
        print(f"Saved to {dirname}/{match}.csv")
    with open(f"{dirname}/{match}_fit.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["a[i]","error[i]"])
        writer.writerows(zip(a,error))
        print(f"Saved to {dirname}/{match}_fit.csv")
    
def extrapolation_plot(match, aw_parameters):
    dirname = aw_parameters["dirname"]
    data = np.genfromtxt(f"{dirname}/{match}.csv",delimiter=',')[1:]
    fig,ax = plot_create()
    ax.plot(data[:,1], data[:,2], label=label_for_plot(match)+" (simulated)", color=color_for_plot(match), marker="o", markersize=12, linestyle="none", zorder=2)
    ys_fitted = np.linspace(data[:,1][0],0,100)
    results_fit = np.genfromtxt(f"{dirname}/{match}_fit.csv",delimiter=',')[1:]
    a = results_fit[:,0]
    error = results_fit[:,1]
    points_fitted = a[0] + a[1]*ys_fitted + a[2]*ys_fitted**2
    ax.plot(ys_fitted, points_fitted, label=label_for_plot(match)+" (fit)", color="black", linewidth=3, zorder=3)
    ax.errorbar(0, a[0], yerr=error[0], color="black", linestyle="", marker="s", capsize=3, barsabove=True, label = "Confidence interval")
    ax.set_xlabel("y = kd|m|")
    ax.legend()
    plot_setaspect(ax)
    plt.savefig(f"{dirname}/{match}.svg", bbox_inches='tight')
    print_log(f"Saved {dirname}/{match}.svg")

def spectrum_with_extrapolation_execute(aw_parameters,adda_cmdlineargs):
    dirname = aw_parameters["dirname"]
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname,exist_ok=True)
    os.makedirs(f"{dirname}/ADDA_output", exist_ok=True)
    print()
    print_log("--- Spectrum with extrapolation: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    mp_files = aw_parameters["mp_files"]
    evs = [x for x in aw_parameters["ev_range"]]
    mdata = dict()
    for i in range(len(mp_files)):
        mdata[i] = np.genfromtxt(mp_files[i],delimiter=',')
    mps = dict()
    m_abs = 0
    for ev in evs:
        mps[ev] = ""
        for i in mdata:
            nk = min(mdata[i], key=lambda x: abs(x[0] - ev))
            mre = float(nk[1])
            mim = float(nk[2])
            if m_abs < math.sqrt(mre**2 + mim**2):
                m_abs = math.sqrt(mre**2 + mim**2)
                maxmre = mre
                maxmim = mim
            nk = str(nk[1]) + " " + str(nk[2]) + " "
            mps[ev] += str(nk)
    size = adda_cmdlineargs["size"]
    grid = adda_cmdlineargs["grid"]
    ev = evs[0]
    lam = ev_to_nm(ev)
    y_min = (2*math.pi/lam)*(size/grid)*m_abs #y = k*d*|m|
    y_max = 4*y_min
    ys = np.exp(np.linspace(math.log(y_min), math.log(y_max), 9))
    grids = np.rint((2*math.pi/lam)*(size/ys)*m_abs).astype(np.int32) #using same grids for all ev - they would not change
    #grids = np.flip([220,186,156,132,110,92,78,66,56])
    adda_cmdlineargs["lambda"] = lam
    adda_cmdlineargs["m"] = f"{maxmre} {maxmim}"
    with open(f"{dirname}/ADDA_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    del adda_cmdlineargs["lambda"]
    del adda_cmdlineargs["m"]
    del adda_cmdlineargs["grid"]
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_files: {mp_files}",dirname)
    print_log(f"Varying energy from {evs[0]} to {evs[-1]} eV",dirname)
    print_log(f"Varying grids: {grids}",dirname)
    cmdlines = []
    for grid_i in grids:
        os.mkdir(f"{dirname}/ADDA_output/{grid_i}")
        for ev in evs:
            cmdline_i = cmdline
            dir_fixed = os.path.abspath(dirname + f"/ADDA_output/{grid_i}/{ev}")
            cmdline_i += f' -dir "{dir_fixed}"'
            cmdline_i += f" -grid {grid_i}"
            cmdline_i += " -lambda %s" % ev_to_nm(ev)
            cmdline_i += " -m " + f" {mps[ev]} "
            cmdline_i += " > " + os.devnull
            cmdlines.append(cmdline_i)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)

def spectrum_with_extrapolation_collect(match, aw_parameters):
    dirname = aw_parameters["dirname"]
    with open(f"{dirname}/ADDA_cmdlineargs.csv") as file:
        reader = csv.reader(file)
        adda_cmdlineargs = dict(reader)
    mre = float(adda_cmdlineargs["m"].split(" ")[0])
    mim = float(adda_cmdlineargs["m"].split(" ")[1])
    m_abs = math.sqrt(mre**2 + mim**2)
    grids = np.array(sorted([int(d.name) for d in os.scandir(f"{dirname}/ADDA_output") if d.is_dir()]))
    evs = sorted([float(d.name) for d in os.scandir(f"{dirname}/ADDA_output/{grids[0]}") if d.is_dir()])
    ys = (2*math.pi/float(adda_cmdlineargs["lambda"]))*(float(adda_cmdlineargs["size"])/grids)*m_abs #y = k*d*|m|
    weights = ys**-3
    fit_values = []
    fit_errors = []
    for ev_i in evs:
        values = []
        for grid_j in grids:
            values.append(parse_value(f"{dirname}/ADDA_output/{grid_j}/{ev_i}/CrossSec-Y",match))
        fit,cov = np.polyfit(ys, values, 2, w=weights, cov=True)
        a = np.flip(fit)
        error = 2*np.sqrt(np.flip(np.diag(cov)))
        fit_values.append(a[0])
        fit_errors.append(error[0])
    with open(f"{dirname}/{match}_fit.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["ev",match,"error"])
        writer.writerows(zip(evs,fit_values,fit_errors))
        print(f"Saved to {dirname}/{match}_fit.csv")

def spectrum_with_extrapolation_plot(match,aw_parameters):
    dirname = aw_parameters["dirname"]
    data = np.genfromtxt(f"{dirname}/{match}_fit.csv",delimiter=',')[1:]
    fig,ax = plot_create()
    ax.plot(data[:,0], data[:,1], label=label_for_plot(match), color=color_for_plot(match), linewidth=3)
    ax.fill_between(data[:,0], data[:,1]-data[:,2], data[:,1]+data[:,2], label="Confidence interval", color="blue", alpha=0.2)
    ax.set_xlim([min(data[:,0]),max(data[:,0])])
    ax.legend()
    plot_setaspect(ax)
    plt.savefig(f"{dirname}/{match}_fit.svg", bbox_inches='tight')
    print_log(f"Saved {dirname}/{match}_fit.svg")

def scan_execute(aw_parameters,adda_cmdlineargs):
    aw_parameters, adda_cmdlineargs = dict(aw_parameters), dict(adda_cmdlineargs)
    dirname = aw_parameters["dirname"]
    start_time = time.time()
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname, exist_ok=True)
    os.makedirs(f"{dirname}/ADDA_output", exist_ok=True)
    print()
    print_log("--- Scan: executing simulations",dirname)
    print_log(f"Number of parallel threads: {aw_parameters['parallel_procs']}",dirname)
    mp_files = aw_parameters["mp_files"]
    if "ev" in aw_parameters.keys():
        ev = aw_parameters["ev"]
        adda_cmdlineargs["lambda"] = ev_to_nm(ev)
    else:
        ev = aw_parameters["nm"]
        adda_cmdlineargs["lambda"] = ev
    mdata = dict()
    for i in range(len(mp_files)):
        mdata[i] = np.genfromtxt(mp_files[i],delimiter=',')
    mps = dict()
    mps[ev] = ""
    for i in mdata:
        nk = min(mdata[i], key=lambda x: abs(x[0] - ev))
        nk = str(nk[1]) + " " + str(nk[2]) + " "
        mps[ev] += str(nk) + " "
    adda_cmdlineargs["m"] = f"{mps[ev]}"
    if 'ms_files' in aw_parameters.keys():
        surfheight = adda_cmdlineargs["surf"].split(" ")[0]
        delkeys_silent(adda_cmdlineargs, ["surf"])
        msdata = np.genfromtxt(aw_parameters["ms_files"][0],delimiter=',')
        nk = min(msdata, key=lambda x: abs(x[0] - ev))
        nk = str(nk[1]) + " " + str(nk[2])
        adda_cmdlineargs["surf"] = f"{surfheight} {nk}"
    size = adda_cmdlineargs["size"]
    grid = adda_cmdlineargs["grid"]
    x_left,x_right = aw_parameters["scan_x_range"]
    y_bottom,y_top = aw_parameters["scan_y_range"]
    step = aw_parameters["scan_step"]
    #adjusting area so the points are exactly in the middle between the dipoles
    d = size/grid #nm
    print_log(f"dipole size = {d} nm",dirname)
    odd = 0.5*(grid % 2)
    x0s, y0s = dipole_middles((x_left,y_bottom),(x_right,y_top),d,odd)
    x0s, y0s = x0s[0::step], y0s[0::step]
    # print(x0s)
    # print(y0s)
    beam_list = adda_cmdlineargs["beam_center"].split(" ")
    del adda_cmdlineargs["beam_center"]
    cmdline = cmdline_construct(aw_parameters,adda_cmdlineargs)
    with open(f"{dirname}/ADDA_cmdlineargs.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        for key, value in adda_cmdlineargs.items():
           writer.writerow([key, value])
    print_log(f"{cmdline}",dirname)
    print_log(f"mp_files: {mp_files}",dirname)
    if "ev" in aw_parameters.keys():
        print_log(f"ev = {ev}",dirname)
    else:
        print_log(f"nm = {ev}",dirname)
    
    print_log(f"m_p = {mps[ev]}",dirname)
    print_log(f"Varying (x_left,x_right) = ({x_left},{x_right}) nm",dirname)
    print_log(f"Varying (y_bottom,y_top) = ({y_bottom},{y_top}) nm",dirname)
    cmdlines = []
    for x0_i in x0s:
        for y0_i in y0s:
            cmdline_i = cmdline
            dir_fixed = os.path.abspath(dirname + f"/ADDA_output/{x0_i}_{y0_i}")
            cmdline_i += f' -dir "{dir_fixed}"'
            beam_list[0], beam_list[1] = str(x0_i), str(y0_i)
            beam_center = (" ").join(beam_list)
            cmdline_i += f" -beam_center {beam_center}"
            cmdline_i += " > " + os.devnull
            cmdlines.append(cmdline_i)
    exec_cmdlines(cmdlines,aw_parameters["parallel_procs"])
    print_log("--- %s seconds" % round((time.time() - start_time),2),dirname)

def scan_collect(match, aw_parameters):
    aw_parameters = dict(aw_parameters)
    dirname = aw_parameters["dirname"]
    dirs = sorted([d.name for d in os.scandir(f"{dirname}/ADDA_output") if d.is_dir()])
    xs = []
    ys = []
    values = []
    for dir_i in dirs:
        xy = dir_i.split("_")
        xs.append(float(xy[0]))
        ys.append(float(xy[1]))
        values.append(parse_value(f"{dirname}/ADDA_output/{dir_i}/CrossSec-Y",match))
    with open(f"{dirname}/{match}.csv", 'w') as file:
        writer = csv.writer(file, delimiter=',', lineterminator='\n')
        writer.writerow(["x","y",match])
        writer.writerows(zip(xs,ys,values))
        print(f"Saved to {dirname}/{match}.csv")

def scan_plot(match, aw_parameters, details=True):
    aw_parameters = dict(aw_parameters)
    dirname = aw_parameters["dirname"]
    with open(f"{dirname}/ADDA_cmdlineargs.csv") as file:
        reader = csv.reader(file)
        adda_cmdlineargs = dict(reader)
    size = float(adda_cmdlineargs["size"])
    grid = float(adda_cmdlineargs["grid"])
    data = np.genfromtxt(f"{dirname}/{match}.csv",delimiter=',',dtype=None, encoding=None)
    axNames = data[0].astype("str")
    data = data[1:].astype("float")
    #print(axNames)
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
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plot_setrcparams()
    ax.set_aspect('equal')
    d = size/grid
    from matplotlib.colors import LogNorm
    im = ax.imshow(z, extent=(min(xs)-d/2, max(xs)+d/2, min(ys)-d/2, max(ys)+d/2), origin="lower", cmap="rainbow")#, norm=LogNorm(vmin=0.01, vmax=1))
    #plt.scatter(x, y, c=z, marker="s") # scatter is the most stable function for visualization, so use this for debugging
    ax.set_xlabel(f"{axNames[0]}, nm")
    ax.set_ylabel(f"{axNames[1]}, nm")
    
    if details == True:
        cbar = fig.colorbar(im)
        cbar.set_label(label_for_plot(axNames[2]))
        #cbar.formatter.set_powerlimits((0, 0))
    else:
        plt.axis('off')

    plt.savefig(f"{dirname}/{match}.png", bbox_inches='tight',dpi=600)
    print_log(f"Saved {dirname}/{match}.png")

