# This code uses amplitude matrix at forward direction to obtain extinction cross sections for various incident
# polarizations and compares them with direct calculations. Also produces circular polarizations from the linear ones
# (relevant for circular dichroism).

import os, shutil, re, numpy as np

fdiff = 1e-8 # fixed relative difference (should be smaller than '-eps 10' used in ADDA runs)

# path to adda executable
#adda_exec = "../../win64/adda.exe"
adda_exec = os.path.abspath(__file__ + "/../../../src/seq/adda")
dir_name = 'out'
beam_name = 'tmpIncBeam'    

def dname(mode):
    return dir_name + str(mode)
    
def str_to_cmplx_list(str_list):
    # builds complex list from string data of floats
    assert len(str_list) % 2 == 0, "Input list must have an even length."
    return np.array([complex(float(str_list[i]), float(str_list[i+1])) for i in range(0, len(str_list), 2)])
    
def cmplx_to_str_list_norm2(cmplx_list):
    # transforms complex list to strings with a squared norm before
    return [str(np.linalg.norm(cmplx_list)**2)] + [str(part) for z in cmplx_list for part in [np.real(z),np.imag(z)]]

def extract_CrossSec(mode,pol):
    # extracts Cext from output
    with open(dname(mode) + '/CrossSec-' + pol) as f:
        file = f.read()
    f.close()
    dat = re.split('[\n \t]',file)
    Cext = float(dat[2])
    return Cext

def extract_ampl(mode):
    # extracts amplitude matrix at forward scattering from output 
    with open(dname(mode) + '/ampl') as f:
        file = f.read()
    f.close()
    dat = re.split('[\n \t]',file)
    return str_to_cmplx_list(dat[10:18])

def adda_run(mode,option):                                                            
    # data generation (run of ADDA code)
    cmdline = ' '.join((adda_exec,'-dir',dname(mode),option,'>',os.devnull))
    os.system(cmdline)

def print_diff(val,x,y):
    # print difference
    cdiff = np.fabs((x-y)/y) # calculated relative difference
    if (cdiff > fdiff):
        print('\n\t\t'+val+':\n\t\tcase 1:\t'+str(x)+'\n\t\tcase 2:\t'+str(y)+'\n\t\tdiff:\t'+str(cdiff)+'\n')
        return 1
    return 0

print('Compare direct calculation of Cext for linear and circular polarization with that through the amplitude matrix.')
# first, standard run with linear polarizations
shape = '-size 4 -shape read helix.dat -eps 10'
adda_run(1,' '.join((shape,'-ntheta 1','-scat_matr ampl','-store_beam')))
CextX = extract_CrossSec(1,'X')
CextY = extract_CrossSec(1,'Y')
S1,S2,S3,S4 = extract_ampl(1)
# in the following k=1 is assumed (as in ADDA runs)
d = 0
d += print_diff('CextX',CextX,4*np.pi*np.real(S1))
d += print_diff('CextY',CextY,4*np.pi*np.real(S2))

# produce left and right circular polarizations from two linear ones, saved during the previous simulation
beamR = beam_name + '-R'
beamL = beam_name + '-L'
with open(dname(1) + '/IncBeam-X', 'r') as fileX, open(dname(1) + '/IncBeam-Y', 'r') as fileY, open(beamR, 'w') as fileR, open(beamL, 'w') as fileL:
    first_line = True
    # Iterate through lines in both input files simultaneously
    for lineX, lineY in zip(fileX, fileY):
        if (first_line):
            # just copy the header line
            assert lineX == lineY, "Header lines for two polarizations must be the same."
            fileR.write(lineX)
            fileL.write(lineY)
            first_line = False
        else:
            datX = lineX.split()
            fieldX = str_to_cmplx_list(datX[4:10])
            datY = lineY.split()
            fieldY = str_to_cmplx_list(datY[4:10])
            # test that the coordinates agree
            assert datX[0:3] == datY[0:3], "Dipole coordinates for two polarizations must be the same."
            # the following definitions are according to Bohren & Huffman, considering that Y and X are parallel and perpendicular polarizations
            # the handiness of the circular polarization is then considered as observed from the detector
            fieldR = (fieldY + 1j*fieldX)/np.sqrt(2)
            fieldL = (fieldY - 1j*fieldX)/np.sqrt(2)
            # The norm of the fields is not exactly 1 due to the limited precision of the original data
            fileR.write(' '.join(datX[0:3] + cmplx_to_str_list_norm2(fieldR)) + '\n')
            fileL.write(' '.join(datX[0:3] + cmplx_to_str_list_norm2(fieldL)) + '\n')
fileX.close()
fileY.close()
fileR.close()
fileL.close()

# second run with circular polarizations, Y and X corresponds to L and R
# iteration method is changed, since circular polarizations seem to cause problems for QMR_CS 
adda_run(2,' '.join((shape,'-scat_matr none','-beam read',beamL,beamR,'-iter bicgstab')))
CextR = extract_CrossSec(2,'X')
CextL = extract_CrossSec(2,'Y')
d += print_diff('CextR',CextR,2*np.pi*(np.real(S1+S2)+np.imag(S4-S3)))
d += print_diff('CextL',CextL,2*np.pi*(np.real(S1+S2)-np.imag(S4-S3)))

# final test resutl
if (d == 0):
    print('\033[32m', '\nPassed', '\033[0m', sep='') #green stdout
else:
    print('\033[31m', '\nNot passed', '\033[0m', sep='') #red

# cleaning
shutil.rmtree(dname(1))
shutil.rmtree(dname(2))
os.remove(beamL)
os.remove(beamR)
