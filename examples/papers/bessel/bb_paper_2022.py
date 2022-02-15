'''
# This code represents the comparison of scattering intensities calculated with GLMT (reference data) and DDA (with ADDA)
# for the scattering of Bessel beam (CS type) by a sphere (option 0) and coated sphere (option 1).
# Reference data were kindly provided by Zhuyang Chen (DOI:10.1088/2040-8978/16/5/055701).
'''

import os, re, math
import matplotlib.pyplot as plt
from matplotlib import rc
from PIL import Image


# path to adda executable
#adda_exec = "../../win64/adda.exe"
adda_exec = os.path.abspath(__file__ + "/../../../../src/seq/adda")


option = 1 # 0 - scattering by a sphere; 1 - scattering by a coated sphere



# Bessel beam parameters
angle = 15 # half-cone angle
lmbd = 0.6328 # beam wavelenght

# particle parameters
eq_rad = lmbd # sphere diameter
grid = 32 # number of dipoles per particle length (see ADDA manual)

run = 1 # 1 - run ADDA code; other - do not run ADDA code


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Data generation (run of ADDA code)
def adda_run(mode): # option 0 or 1
    
    # Beam and particle definition

    run_options = [' -beam besselCS 0 ' + str(angle) + ' -shape coated 0.5 -m 1.33 0 1.33 0', # option 0 - sphere
                   ' -beam besselCS 0 ' + str(angle) + ' -shape coated 0.5 -m 1.33 0 1.55 0'  # option 1 - coated sphere
                   ]                                                               
    
    # cmd line generation (see ADDA manual)
    cmdline = adda_exec
    cmdline += ' -sym enf' #Do not simulate second polarization
    cmdline += ' -ntheta 180'
    cmdline += run_options[mode]
    cmdline += ' -lambda ' + str(lmbd)
    cmdline += ' -store_beam'
    cmdline += ' -eq_rad ' + str(eq_rad)
    cmdline += ' -grid ' + str(grid)
    cmdline += ' -dir dda/option_' + str(mode)
    print(cmdline)
    
    os.system(cmdline)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data extraction
def extractData(mode):
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
    
    path = 'dda/option_' + str(mode)
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
    plt.plot(xv1, yv1, label = 'DDA', color = (0.5, 0.7, 0.9))
    plt.plot(xv2, yv2, label = 'GLMT', color = 'm',linestyle='dashed')
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


# images of the scattering particle shape
def plotShape(size):
    ax.set_title(r'D = ' + str(size) + ' $\mu m$')
    img = Image.open('dda/option_'+str(option)+'/shape_image.png')
    ax.set_axis_off()
    ax.imshow(img)


# Visualisation of the amplitude of the incident electric field almost in the middle
# of the particle (z = distance to the closest dipole along z axis)
def plotField(xd,yd,ed,z0):
    ax.set_title('Amplitude of the incident field \n(z = '+str(round(z0,2))+')');
    ax.scatter(xd, yd, ed, c=ed, cmap='viridis', linewidth=0.5);
    #ax.plot_trisurf(xd, yd, ed,cmap='viridis', edgecolor='none');
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# whether run ADDA code or not
if run == 1:
    adda_run(option)
    
dda_theta,dda_iper,dda_ipar,dda_xd,dda_yd,dda_ed,dda_z0 = extractData(option) # extraction of ADDA results
ref_thetaper,ref_iper,ref_thetapar,ref_ipar = extractData('glmt') # extraction of reference data (GLMT)


# data visualisation

fig = plt.figure()
# Image of particle shape
ax = fig.add_subplot(221)
plotShape(eq_rad*2)
# Visualisation of the amplitude of the incident electric field
ax = fig.add_subplot(222,projection='3d')
plotField(dda_xd,dda_yd,dda_ed,dda_z0)
# Perpendicular scattering intensity (1)
ax = fig.add_subplot(223)
plotData(dda_theta,dda_iper,ref_thetaper,ref_iper,1)
# Parallel scattering intensity (2)
ax = fig.add_subplot(224)
plotData(dda_theta,dda_ipar,ref_thetapar,ref_ipar,2)

plt.tight_layout()

# data save
os.makedirs('saved', exist_ok=True)
if option == 0:
    plt.savefig('saved/BB_sphere.pdf', bbox_inches='tight')
if option == 1:
    plt.savefig('saved/BB_coatedsphere.pdf', bbox_inches='tight')
    
plt.show()