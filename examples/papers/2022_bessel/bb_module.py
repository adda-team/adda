'''
# The module contains functions for the demonstration of figures 12-15 in the 
# corresponding python scripts (fig... .py). 
'''

import os, shutil, re, math
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np


# path to adda executable
#adda_exec = "../../win64/adda.exe"
adda_exec = os.path.abspath(__file__ + "/../../../../src/seq/adda")


fig1214comm = ' -beam besselCS 0 15 -eq_rad 0.6328 -lambda 0.6328'
fig15comm = ' -m 1.52 0 -shape box -size 1. -store_beam -lambda 0.6328 -grid 15'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Data generation (run of ADDA code)
def adda_run(mode):
    
    if os.path.exists('dda/option_' + str(mode)):
        shutil.rmtree('dda/option_' + str(mode))
    
    # Beam and particle definition
    run_options = [fig1214comm + ' -m 1.33 0 -grid 64', # option 0 - sphere
                   fig1214comm + ' -shape coated 0.5 -m 1.33 0 1.55 0 -grid 64',  # option 1 - coated sphere
                   fig15comm + ' -beam besselCS  4 45',
                   fig15comm + ' -beam besselCSp 4 45',
                   fig15comm + ' -beam besselTEL 4 45',
                   fig15comm + ' -beam besselTML 4 45'
                   ]                                                               
    
    # cmd line generation (see ADDA manual)
    cmdline = adda_exec
    cmdline += ' -sym enf' #Do not simulate second polarization
    cmdline += ' -ntheta 180'
    cmdline += run_options[mode]
    cmdline += ' -dir dda/option_' + str(mode)
    print(cmdline)
    
    os.system(cmdline + ' > ' + os.devnull)
    
    
def adda_run_grids(grids): # mode from 0 up to 5      
    # Beam and particle definition
    run_option = fig1214comm + ' -m 1.33 0'
    # cmd line generation
    cmdline = adda_exec
    cmdline += ' -sym enf' #Do not simulate second polarization
    cmdline += ' -ntheta 90'
    cmdline += run_option
    
    for g in grids:
        if os.path.exists('dda/extrapolation/run_' + str(g)):
            shutil.rmtree('dda/extrapolation/run_' + str(g))
        os.system(cmdline + ' -grid ' + str(g) + ' -dir dda/extrapolation/run_' + str(g))
    return 0


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data extraction
def extractData(mode,option):
    if mode == 'glmt':
        # Reference data (Generalized Lorenz Mie Theory)
        ref_path = 'ref/glmt/option_' + str(option) + '/'
        files_ref = os.listdir(ref_path)
        # Scattering in perpendicular plane
        with open(ref_path + str(files_ref[files_ref.index('glmt1.txt')])) as f:
            file_ref = f.read()
        f.close()
        dat_ref = re.split('[\n \t]',file_ref)
        theta_per = dat_ref[0:-1:2]
        theta_per = [float(i) for i in theta_per]
        i_per = dat_ref[1::2]
        i_per = [math.pi*pow(10,float(i)/10) for i in i_per]
        # Scattering in parallel plane
        with open(ref_path + str(files_ref[files_ref.index('glmt2.txt')])) as f:
            file_ref = f.read()
        f.close()
        dat_ref = re.split('[\n \t]',file_ref)
        theta_par = dat_ref[0:-1:2]
        theta_par = [float(i) for i in theta_par]
        i_par = dat_ref[1::2]
        i_par = [math.pi*pow(10,float(i)/10) for i in i_par]
        
        return theta_per,i_per,theta_par,i_par
    
    path = 'dda/option_' + str(option)
    files = os.listdir(path)
    with open(path + '/' + str(files[files.index('mueller')])) as f:
        fileM = f.read()
    f.close()
    dat = re.split('[\n \t]',fileM)
    theta = dat[17:-1:17]
    s11 = dat[18:-1:17]
    s12 = dat[19:-1:17]
    
    theta = [float(th) for th in theta]
    s11 = [float(s) for s in s11]
    s12 = [float(s) for s in s12]
    i_per = [x-y for x, y in zip(s11, s12)]
    i_par = [x+y for x, y in zip(s11, s12)]
    
    return theta,i_per,i_par

# fig 15 intensity profiles
def extractField(mode):
    path = 'dda/option_' + str(mode)
    files = os.listdir(path)
    with open(path + '/' + str(files[files.index('IncBeam-Y')])) as f:
        fileF = f.read()
    f.close()
    dat = re.split('[\n \t]',fileF)
    xd = dat[10:-1:10]
    yd = dat[11:-1:10]
    zd = dat[12:-1:10]
    ed = dat[13:-1:10]
    xd = [float(i) for i in xd]
    yd = [float(i) for i in yd]
    zd = [float(i) for i in zd]
    ed = [float(i) for i in ed]
    minz = min([abs(i) for i in zd])
    zslice = [i for i in range(len(zd)) if zd[i]==minz]
    xd = [xd[i] for i in zslice]
    yd = [yd[i] for i in zslice]
    ed = [ed[i] for i in zslice]
    
    return xd,yd,ed,minz

# fig 14 extrapolation
def extractData_grids(grids,tht):
    iperg = []
    for g in grids:
        path = 'dda/extrapolation/run_' + str(g)
        files = os.listdir(path)
        with open(path + '/' + str(files[files.index('mueller')])) as f:
            fileM = f.read()
        f.close()
        dat = re.split('[\n \t]',fileM)
        thts = dat[17:-1:17]
        s11 = dat[18:-1:17]
        s12 = dat[19:-1:17]
        thts = [float(th) for th in thts]
        s11 = [float(s) for s in s11]
        s12 = [float(s) for s in s12]
        i_per = [x-y for x, y in zip(s11, s12)]
        ind = thts.index(tht)
        iperg.append(i_per[ind])
    
    return iperg

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data visualisation

SMALL_SIZE = 10
MEDIUM_SIZE = 12

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font', size=MEDIUM_SIZE) # controls default text sizes

# fig 12-13 scattering intensities
def plotData(xv1,yv1,xv2,yv2,flag,ax):
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    plt.plot(xv1, yv1, label = 'DDA', color = (107/256,18/256,245/256),linewidth = 1)
    plt.plot(xv2, yv2, label = 'GLMT', color = (234/256,56/256,38/256),linestyle='none',marker = '.',markersize = 3)
    plt.minorticks_on()
    plt.tick_params(which='major',right=True, top=True)
    plt.tick_params(which='minor',right=True, top=True)
    plt.xlabel(r'Scattering angle $\theta$ [deg]',labelpad=5.)
    if flag == 1:
        plt.ylabel(r'Parallel intensity $I_{\parallel}$',labelpad=2.)
        plt.title('(a)',loc='left', x=0.035,y=0.03)
    if flag == 2:
        plt.ylabel(r'Perpendicular intensity $I_{\perp}$',labelpad=2.)
        plt.title('(b)',loc='left', x=0.035,y=0.03)
    plt.legend()
    plt.yscale('log')
    plt.xlim(0, 180)
    plt.xticks(np.arange(0, 181, 30))
    x_left1, x_right1 = ax.get_xlim()
    y_low1, y_high1 = ax.get_ylim()
    ax.set_aspect(abs((x_right1-x_left1)/(math.log(y_low1)-math.log(y_high1)))*1.85,adjustable="box")

# fig 15 scattering intensities
def plotData4(xv1,yv1,xv2,yv2,xv3,yv3,xv4,yv4,flag,ax):
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    plt.plot(xv1, yv1, label = r'CS', color = (107/256,18/256,245/256), linewidth=1.5,linestyle='dashed')
    plt.plot(xv2, yv2, label = "CS'", color = (121/256,251/256,186/256),linewidth=1.5,linestyle='solid')
    plt.plot(xv3, yv3, label = 'TEL', color = (230/256,139/256,79/256), linewidth=1.5,linestyle='dotted')
    plt.plot(xv4, yv4, label = 'TML', color = (234/256,56/256,38/256), linewidth=1.5,linestyle='dashdot')
    plt.minorticks_on()
    plt.tick_params(which='major',right=True, top=True)
    plt.tick_params(which='minor',right=True, top=True)
    plt.xlabel(r'Scattering angle $\theta$ [deg]',labelpad=5.)
    if flag == 1:
        plt.ylabel(r'Parallel intensity $I_{\parallel}$',labelpad=3.)
        plt.title('(e)',loc='left', x=0.035,y=0.85)
    if flag == 2:
        plt.ylabel(r'Perpendicular intensity $I_{\perp}$',labelpad=3.)
        plt.legend(loc='lower left',ncol=2,fontsize=SMALL_SIZE)
        plt.title('(f)',loc='left', x=0.035,y=0.85)
    plt.yscale('log')
    plt.xlim(0,180)
    plt.xticks(np.arange(0, 181, 30))
    x_left1, x_right1 = ax.get_xlim()
    y_low1, y_high1 = ax.get_ylim()
    ax.set_aspect(abs((x_right1-x_left1)/(math.log(y_low1)-math.log(y_high1)))*1.65,adjustable="box")
    
# fig 15
# Visualisation of the intensity of the incident electric field almost in the 
# middle of the particle (the nearest to the center xy-plane of dipoles is 
# used, its z-coordinate is shown on the plot)
def plotField(xd,yd,ed,z0,mode,ax):
    axtitles = ['(a)',"(b)",'(c)','(d)']
    axtypes = ['CS',"CS'",'TEL','TML']
    ax.set_title(axtitles[mode-2],loc='left',y=0.8,x=0.08)
    ax.scatter(xd, yd, ed, c=ed, cmap="rainbow",label=axtypes[mode-2])
    plt.legend(loc=[0.1,0.1],frameon=False,markerfirst=False,markerscale=0)
    ax.axis('off')

# fig 14 extrapolation
def axesGen(grd):
    #grids_ut = grd[N0:N1:]
    grids_ut = grd
    lmbd = 0.6328
    m_abs = 1.33
    size = lmbd*2
    ys = list(map(lambda x: (2*math.pi/lmbd)*(size/x)*m_abs, grids_ut))
    #ys = list(map(lambda x: (1/x), grids_ut))
    return ys

# fig 14 extrapolation
def fitData(grids,data):
    #values = data[N0:N1:]
    values = data
    ys = axesGen(grids)
    weights = list(map(lambda x: x**-3, ys))
    fit,cov = np.polyfit(ys, values, 2, w=weights, cov=True)
    a = np.flip(fit)
    error = 2*np.sqrt(np.flip(np.diag(cov)))
    return a,error

# fig 14 extrapolation
def fitDotPlot(grids,data,ref):
    a,error = fitData(grids,data)
    ys = axesGen(grids)
    ymax = 1/(1/(ys[0]+0.1))
    y = np.arange(0,ymax,ymax/200)
    res_ext = list(map(lambda x: a[0]+x*a[1] + (x**2)*a[2], y))

    plt.plot(y, res_ext, label = 'Extrapolation', color = 'gray')
    plt.plot(ys, data, label = 'DDA', marker="o", linestyle="none",markersize = 6, color = (124/256,5/256,254/256))
    plt.errorbar(0, a[0], yerr=error[0], color="k", marker="s", capsize=1, barsabove=True)
    plt.plot([0], ref, label = 'GLMT', marker="d", markersize = 6, linestyle="none", color = (255/256,22/256,11/256))
    plt.minorticks_on()
    plt.ylabel('Backscattering intensity')
    plt.xlabel(r'Discretization parameter $kdn_1$')
    plt.tick_params(which='major',right=True, top=True)
    plt.tick_params(which='minor',right=True, top=True)
    plt.legend(loc='lower right')
    return 0