'''
# This code compares the scattering intensities calculated with the DDA (ADDA code)
# for the scattering of two different Bessel beams (option_1 and option_2) by a sphere.
'''

import os, re, math
import matplotlib.pyplot as plt
from matplotlib import rc


# path to adda executable
#adda_exec = "../../win64/adda.exe"
adda_exec = os.path.abspath(__file__ + "/../../../src/seq/adda")


# define here different parameters for 2 options (see ADDA manual)
run_options = [' -beam besselLE  2 15',  # option 1
               ' -beam besselLM  2 15']  # option 2

# =============================================================================
# Bessel beams in ADDA:
#     
#   besselCS <order> <angle>
#       -Bessel beam with circularly symmetric energy density.
#   besselCSp <order> <angle>
#       -Alternative Bessel beam with circularly symmetric energy density.
#   besselM <order> <angle> <ReMex> <ReMey> <ReMmx> <ReMmy> [<ImMex> <ImMey> <ImMmx> <ImMmy>]
#       -Generalized Bessel beam. The beam is defined by 2x2 matrix M:
#       (Mex, Mey, Mmx, Mmy). Real parts of these four elements are obligatory, 
#       while imaginary parts are optional (zero, by default).
#   besselLE <order> <angle>
#       -Bessel beam with linearly polarized electric field.
#   besselLM <order> <angle>
#       -Bessel beam with linearly polarized magnetic field.
#   besselTEL <order> <angle>
#       -Linear component of the TE Bessel beam.
#   besselTML <order> <angle>
#       -Linear component of the TM Bessel beam.
#   
#   Order is integer (of any sign) and the half-cone angle (in degrees) 
#   is measured from the z-axis.
# =============================================================================
                                                                              

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data generation (run of ADDA code)
def adda_run(mode): # option 1 or 2                                                             
    
    # cmd line generation (see ADDA manual)
    
    # common parameters for 2 options
    cmdline = adda_exec
    cmdline += ' -grid 16' # particle discretization
    cmdline += ' -sym enf' # do not simulate second polarization
    cmdline += ' -ntheta 180' # number of scattering angles
    cmdline += ' -store_beam' # save incident field
    cmdline += ' -dir option_' + str(mode) # save path
    cmdline += run_options[mode-1]
    print(cmdline)
    
    os.system(cmdline)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data extraction
def extractData(mode):
    path = 'option_' + str(mode)
    files = os.listdir(path)
    with open(path + '/' + str(files[files.index('mueller')])) as f:
        fileM = f.read()
    f.close()
    dat = re.split('[\n \t]',fileM)
    theta = dat[17:-1:17]
    s11 = dat[18:-1:17]
    s12 = dat[19:-1:17]
    
    with open(path + '/' + str(files[files.index('IncBeam-Y')])) as f:
        fileF = f.read()
    f.close()
    dat = re.split('[\n \t]',fileF)
    xd = dat[10:-1:10]
    yd = dat[11:-1:10]
    zd = dat[12:-1:10]
    ed = dat[13:-1:10]
    
    theta = [float(th) for th in theta]
    s11 = [float(s) for s in s11]
    s12 = [float(s) for s in s12]
    i_per = [x-y for x, y in zip(s11, s12)]
    i_par = [x+y for x, y in zip(s11, s12)]
    xd = [float(i) for i in xd]
    yd = [float(i) for i in yd]
    zd = [float(i) for i in zd]
    ed = [float(i) for i in ed]
    minz = min([abs(i) for i in zd])
    zslice = [i for i in range(len(zd)) if zd[i]==minz]
    xd = [xd[i] for i in zslice]
    yd = [yd[i] for i in zslice]
    ed = [ed[i] for i in zslice]
    
    return theta,i_per,i_par,xd,yd,ed,minz

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data visualisation

# plots of scattering intensities
def plotData(xv1,yv1,xv2,yv2,flag):
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    plt.plot(xv1, yv1, label = 'option 1', color = (0.5, 0.7, 0.9))
    plt.plot(xv2, yv2, label = 'option 2', color = 'm',linestyle='dashed')
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


# Visualisation of the amplitude of the incident electric field almost in the middle
# of the particle (the nearest to the center xy-plane of dipoles is used, its z-coordinate is shown on the plot)
def plotField(xd,yd,ed,z0,mode):
    ax.set_title(r'OPTION '+str(mode)+':\n '+run_options[mode-1]+'\nIntensity profile of $|E_{inc}|^2$\n(z = '+str(round(z0,2))+')');
    ax.scatter(xd, yd, ed, c=ed, cmap='viridis', linewidth=0.5);
    #ax.plot_trisurf(xd, yd, ed,cmap='viridis', edgecolor='none');
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

adda_run(1)
adda_run(2)
    
theta_1,iper_1,ipar_1,xd_1,yd_1,ed_1,z0_1 = extractData(1) # extraction of ADDA results for option 1
theta_2,iper_2,ipar_2,xd_2,yd_2,ed_2,z0_2 = extractData(2) # extraction of ADDA results for option 2


# data visualisation

fig = plt.figure()
# Visualisation of the amplitude of the incident electric field for option 1
ax = fig.add_subplot(221,projection='3d')
plotField(xd_1,yd_1,ed_1,z0_1,1)
# Visualisation of the amplitude of the incident electric field for option 2
ax = fig.add_subplot(222,projection='3d')
plotField(xd_2,yd_2,ed_2,z0_2,2)
# Perpendicular scattering intensity (1)
ax = fig.add_subplot(223)
plotData(theta_1,iper_1,theta_2,iper_2,1)
# Parallel scattering intensity (2)
ax = fig.add_subplot(224)
plotData(theta_1,ipar_1,theta_2,ipar_2,2)

plt.tight_layout()


# data save
os.makedirs('saved', exist_ok=True)

plt.savefig('saved/results.pdf', bbox_inches='tight')
        
plt.show()