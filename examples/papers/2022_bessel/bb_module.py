'''
# This code represents the comparison of scattering intensities calculated 
# with GLMT (reference data) and DDA (with ADDA) for the scattering of 
# Bessel beam (CS type) by a sphere (option 0) and coated sphere (option 1).
# Reference data were kindly provided by Zhuyang Chen 
# (DOI:10.1088/2040-8978/16/5/055701).
'''

import os, re, math
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np


# path to adda executable
#adda_exec = "../../win64/adda.exe"
adda_exec = os.path.abspath(__file__ + "/../../../../src/seq/adda")


option = 0 # 0 - scattering by a sphere; 1 - scattering by a coated sphere



# Bessel beam parameters
angle = 15 # half-cone angle
lmbd = 0.6328 # beam wavelenght

# particle parameters
eq_rad = lmbd # sphere diameter
grid = 32 # number of dipoles per particle length (see ADDA manual)
# Constants 
mre = 1.33
mim = 0

run = 1 # 1 - run ADDA code; other - do not run ADDA code


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Data generation (run of ADDA code)
def adda_run(mode):
    
    # Beam and particle definition

    run_options = [' -beam besselCS 0 ' + str(angle) + ' -m 1.33 0 -eq_rad ' + str(eq_rad), # option 0 - sphere
                   ' -beam besselCS 0 ' + str(angle) + ' -shape coated 0.5 -m 1.33 0 1.55 0  -eq_rad ' + str(eq_rad),  # option 1 - coated sphere
                   ' -beam besselCS 2 15 -beam_center 0 0 0 -m 1.33 0  -shape box -size ' + str(2*lmbd),
                   ' -beam besselCS 2 15 -beam_center 0 '+str(lmbd)+'   0 -m 1.33 0 -shape box -size ' + str(2*lmbd),
                   ' -beam besselCS 2 15 -beam_center 0 '+str(2*lmbd)+' 0 -m 1.33 0 -shape box -size ' + str(2*lmbd),
                   ]                                                               
    
    # cmd line generation (see ADDA manual)
    cmdline = adda_exec
    #cmdline += ' -sym enf' #Do not simulate second polarization
    cmdline += ' -ntheta 360'
    cmdline += run_options[mode]
    cmdline += ' -lambda ' + str(lmbd)
    #cmdline += ' -eq_rad ' + str(eq_rad)
    cmdline += ' -grid ' + str(grid)
    cmdline += ' -dir dda/option_' + str(mode)
    print(cmdline)
    
    os.system(cmdline + ' > ' + os.devnull)
    
    
def adda_run_grids(grids): # mode from 0 up to 5
    # Beam and particle definition
    run_option = ' -beam besselCS 0 ' + str(angle) + ' -m 1.33 0 -eq_rad ' + str(eq_rad)
    # cmd line generation
    cmdline = adda_exec
    cmdline += ' -sym enf' #Do not simulate second polarization
    cmdline += ' -ntheta 90'
    cmdline += run_option
    cmdline += ' -lambda ' + str(lmbd)
    
    for g in grids:
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

# plots of scattering intensities
def plotData(xv1,yv1,xv2,yv2,flag,ax):
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    plt.plot(xv1, yv1, label = 'DDA', color = (107/256,18/256,245/256))
    plt.plot(xv2, yv2, label = 'GLMT', color = (234/256,56/256,38/256),linestyle='none', marker = '.', markersize = 3)
    plt.minorticks_on()
    plt.tick_params(which='major',right=True, top=True)
    plt.tick_params(which='minor',right=True, top=True)
    plt.xlabel(r'Scattering angle $\theta$, deg')
    if flag == 1:
        plt.ylabel(r'$I_{\perp}$ (H plane)')
    if flag == 2:
        plt.ylabel(r'$I_{\parallel}$ (E plane)')
    plt.legend()
    plt.xlim(0, 180)
    plt.yscale('log')
    x_left1, x_right1 = ax.get_xlim()
    y_low1, y_high1 = ax.get_ylim()
    ax.set_aspect(abs((x_right1-x_left1)/(math.log(y_low1)-math.log(y_high1)))*2,adjustable="box")
    '''
    maxint = max(yv2)
    dev = list(map(lambda x, y: math.fabs((x-y)/maxint*100), yv1,yv2))
    dev2 = list(map(lambda x, y: (math.fabs((x-y)/maxint))**2, yv1,yv2))

    print('\nDeviations ' + str(flag))
    print('Max relative error, %:')
    print(max(dev))
    print('Average relative error,%:')
    print(sum(dev)/len(dev))
    print('RMS, %:')
    print(math.sqrt(sum(dev2)/len(xv2)*100))
    '''

def plotData3(xv1,yv1,xv2,yv2,xv3,yv3,flag,ax):
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    plt.plot(xv1, yv1, label = r'$\Delta$y = 0', color = (107/256,18/256,245/256), linewidth = 1)
    plt.plot(xv2, yv2, label = r'$\Delta$y = a', color = (121/256,251/256,186/256),linewidth=1)
    plt.plot(xv3, yv3, label = r'$\Delta$y = 2a', color = (234/256,56/256,38/256), linewidth = 1)
    plt.minorticks_on()
    plt.tick_params(which='major',right=True, top=True)
    plt.tick_params(which='minor',right=True, top=True)
    plt.xlabel(r'Scattering angle $\theta$, deg')
    if flag == 1:
        plt.ylabel(r'$I_{\perp}$ (H plane)')
    if flag == 2:
        plt.ylabel(r'$I_{\parallel}$ (E plane)')
        plt.legend(loc='upper center')
    plt.xlim(0, 360)
    plt.yscale('log')
    x_left1, x_right1 = ax.get_xlim()
    y_low1, y_high1 = ax.get_ylim()
    ax.set_aspect(abs((x_right1-x_left1)/(math.log(y_low1)-math.log(y_high1)))*2,adjustable="box")
    
    
    
def axesGen(grd):
    #grids_ut = grd[N0:N1:]
    grids_ut = grd
    m_abs = math.sqrt(mre**2 + mim**2)
    size = eq_rad*2
    ys = list(map(lambda x: (2*math.pi/lmbd)*(size/x)*m_abs, grids_ut))
    #ys = list(map(lambda x: (1/x), grids_ut))
    return ys

def fitData(grids,data):
    #values = data[N0:N1:]
    values = data
    ys = axesGen(grids)
    weights = list(map(lambda x: x**-3, ys))
    fit,cov = np.polyfit(ys, values, 2, w=weights, cov=True)
    a = np.flip(fit)
    error = 2*np.sqrt(np.flip(np.diag(cov)))
    return a,error

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