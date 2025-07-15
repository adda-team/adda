#!/usr/bin/env python3

# This code uses amplitude matrix at forward direction to obtain extinction cross sections for various incident
# polarizations and compares them with direct calculations. Also produces circular polarizations from the linear ones
# (relevant for circular dichroism). Finally, tests different options to emulate circular polarizations using
# plane wave limit of Bessel beams.
#
# TODO: Some of the functions are taken from bb_equiv.py and can be reused as module, 
# but probably it is better to take them from ADDAwrapper


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
# in the following k=1 is assumed (as in ADDA runs) - see Section 11.4 "Integral scattering quantities" of the manual
CextX_ampl = 4*np.pi*np.real(S1)
CextY_ampl = 4*np.pi*np.real(S2)
CextR_ampl = 2*np.pi*(np.real(S1+S2)+np.imag(S4-S3))
CextL_ampl = 2*np.pi*(np.real(S1+S2)-np.imag(S4-S3))
# compare Cext for linear polarizations
d = 0
d += print_diff('CextX',CextX,CextX_ampl)
d += print_diff('CextY',CextY,CextY_ampl)

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
            # The following definitions are according to Bohren & Huffman, considering that Y and X are parallel and
            # perpendicular polarizations. The handiness of the circular polarization is then considered as observed
            # from the detector
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
d += print_diff('CextR',CextR,CextR_ampl)
d += print_diff('CextL',CextL,CextL_ampl)

# We can also produce circular polarizations through Bessel beams, first LE^(1,{+,-}i) 
# We use the following from (Glukhova & Yurkin, 2022):
#   x -> y (perpendicular to parallel transformation): Eq.(53)
#   R.e_{+,-} = {-,+}i*e_{+,-}: p.10, second column, first line
#   e_{+,-} = {+,-}i*sqrt(2)*e_{L,R}: p.3, first column, last line
#   limit of LE^(x) is ex: Table I.
# Thus, LE^(1,{+,-}i) goes to e_{+,-}, the corresponding Y-polarization (parallel) is then sqrt(2)*e_{L,R}
# Scaling the M matrix (Tabe II) by 1/sqrt(2), we get exactly e_{L,R} as incident field.
# Here and furher, note that switch from R to L is equivalent to conjugating M. 
line_except_beam = ' '.join((shape,'-scat_matr none','-sym enf','-iter bicgstab'))
tmp = 2**(-1/2)
adda_run(3,' '.join((line_except_beam,'-beam besselM 0 0 0 0 0',str(tmp),'0 0',str(tmp),'0')))
CextR_LE = extract_CrossSec(3,'Y')
adda_run(4,' '.join((line_except_beam,'-beam besselM 0 0 0 0 0',str(tmp),'0 0',str(-tmp),'0')))
CextL_LE = extract_CrossSec(4,'Y')
d += print_diff('CextR_LE',CextR_LE,CextR_ampl)
d += print_diff('CextL_LE',CextL_LE,CextL_ampl)

# Exactly analogous is for LM^(a,b) beams, but the limit of LM^(x) is ey, then LM^(1,{+,-}i) goes to {-,+}i*e_{+,-}.
# So we scale the corresponding M matrix by {+,-}i/sqrt(2) to get e_{L,R}
tmp = 2**(-1/2)
adda_run(5,' '.join((line_except_beam,'-beam besselM 0 0',str(tmp),'0 0 0 0',str(-tmp),'0 0')))
CextR_LM = extract_CrossSec(5,'Y')
adda_run(6,' '.join((line_except_beam,'-beam besselM 0 0',str(tmp),'0 0 0 0',str(tmp),'0 0')))
CextL_LM = extract_CrossSec(6,'Y')
d += print_diff('CextR_LM',CextR_LM,CextR_ampl)
d += print_diff('CextL_LM',CextL_LM,CextL_ampl)

# More intricate approach is through bessel CS' beams with orders 2 and -2. It has shorter command line,
# but needs to be scaled afterwards, since the half-cone angle alpha can not be exactly zero
# We use additionally from (Glukhova & Yurkin, 2022):
#   lim{a->0}E(x)_CS' = a^2*e+/8: p.9, second-to-last line
#   inverting the sign of beam order is equivalent to conjugating everything
tmp = 1e-6 # for scaling half-cone angle, relative errors are expected to be square of that
alpha = (2**(5/4))*(180/np.pi)*tmp
# based on the above, the Y-polarizations of incident beams for the following runs are -tmp^2
# times R and L polarizations, respectively
adda_run(7,' '.join((line_except_beam,'-beam besselCSp -2',str(alpha))))
CextR_CSp = extract_CrossSec(7,'Y')/(tmp**4)
adda_run(8,' '.join((line_except_beam,'-beam besselCSp 2',str(alpha))))
CextL_CSp = extract_CrossSec(8,'Y')/(tmp**4)
d += print_diff('CextR_CSp',CextR_CSp,CextR_ampl)
d += print_diff('CextL_CSp',CextL_CSp,CextL_ampl)

# It can also be done through CS beams (see, again, Table I and II), but the resulting M matrix will be
# average of that for LE and LM beams, so this test is redundant.

# Finally, we do it through TE beams of orders 1 or -1, which are equivalent to TEL^(1,{+,-}i), see Eq.(63).
# We then use the limit in Table I, or Eq.(59), and scale the matrix M from Table II (the benefit is that we do not
# need to scale cross sections afterwards). Half-cone angle should be non-zero, as for CS' above, but the calculation
# inside ADDA is currently not robust to the cancellations, on which we rely here. Thus, the value of tmp below is
# quasi-optimal, since smaller values lead to increasing errors. This value leads to errors around 1e-7, so we need
# to update the tolerance.
# TODO: the following need to be removed, when the code in ADDA is made more robust
fdiff = 1e-6
tmp = 1e-4 # for scaling half-cone angle, relative errors are expected to be square of that
al_rad = np.sqrt(2)*tmp
alpha = al_rad*180/np.pi        # value in degrees for the command line
sc_k = -1/(tmp*np.sin(al_rad))  # scaled k/kt
sc_kz = -1/(tmp*np.tan(al_rad)) # scaled kz/kt
adda_run(9,' '.join((line_except_beam,'-beam besselM 0',str(alpha),str(-sc_k),'0 0',str(sc_kz),'0',str(sc_k),str(sc_kz),'0')))
CextR_TEL = extract_CrossSec(9,'Y')
adda_run(10,' '.join((line_except_beam,'-beam besselM 0',str(alpha),str(-sc_k),'0 0',str(sc_kz),'0',str(-sc_k),str(-sc_kz),'0')))
CextL_TEL = extract_CrossSec(10,'Y')
d += print_diff('CextR_TEL',CextR_TEL,CextR_ampl)
d += print_diff('CextL_TEL',CextL_TEL,CextL_ampl)

# analogously, for TM or TML^(1,{+,-}i). Similarly, as for TE,TM pair above, here we additionally multiply M by {+,-}i for {L,R}
# this also follows from limits in Eq.(59). The result is an interchange of k and kz and change of sign (which is logical, since
# the difference between k and kz determines the limit).
adda_run(11,' '.join((line_except_beam,'-beam besselM 0',str(alpha),str(sc_kz),'0 0',str(-sc_k),'0',str(-sc_kz),str(-sc_k),'0')))
CextR_TML = extract_CrossSec(11,'Y')
adda_run(12,' '.join((line_except_beam,'-beam besselM 0',str(alpha),str(sc_kz),'0 0',str(-sc_k),'0',str(sc_kz),str(sc_k),'0')))
CextL_TML = extract_CrossSec(12,'Y')
d += print_diff('CextR_TML',CextR_TML,CextR_ampl)
d += print_diff('CextL_TML',CextL_TML,CextL_ampl)

# final test of accumulated error statuses
if (d == 0):
    print('\033[32m', '\nPassed', '\033[0m', sep='') #green stdout
else:
    print('\033[31m', '\nNot passed', '\033[0m', sep='') #red

# cleaning
for i in range(1,13): 
    shutil.rmtree(dname(i))
os.remove(beamL)
os.remove(beamR)
