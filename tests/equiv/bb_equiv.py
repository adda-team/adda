'''
# This code compares equivalent cases (command lines)
# for the scattering of Bessel beams calculated in ADDA.
'''

import os,shutil,re,math



fdiff = 1e-10 # fixed difference


# path to adda executable
#adda_exec = "../../win64/adda.exe"
adda_exec = os.path.abspath(__file__ + "/../../../src/seq/adda")
dirname = 'out'    

                                                                    

def extractCrosSeq(mode,pol):
    dname = dirname + str(mode)
    files = os.listdir(dname)
    with open(dname + '/' + str(files[files.index('CrossSec-' + pol)])) as f:
        file = f.read()
    f.close()
    dat = re.split('[\n \t]',file)
    Cext = float(dat[2])
    Qext = float(dat[5])
    Cabs = float(dat[8])
    Qabs = float(dat[11])
    
    return Cext,Qext,Cabs,Qabs


# data generation (run of ADDA code)
def adda_run(mode,option):                                                            
    dname = dirname + str(mode)
    #os.makedirs(dname, exist_ok=True)
    cmdline = adda_exec + ' -dir ' + dname + option + ' > ' + os.devnull
    os.system(cmdline)


# print difference
def printdiff(val,x,y):
    cdiff = math.fabs(x-y) # calculated difference
    if (cdiff > fdiff):
        fr.write('\n\t\t'+val + ':\n\t\tcase 1:\t'+str(x)+'\n\t\tcase 2:\t'+str(y)+'\n\t\tdiff:\t'+str(cdiff)+'\n')
        return 1
    return 0


def printsearch(pol):
    fr.write('\n\n\tCrosSec-'+pol+': search diff')
    cext1,qext1,cabs1,qabs1 = extractCrosSeq(1,pol)
    cext2,qext2,cabs2,qabs2 = extractCrosSeq(2,pol)
    d = 0
    d += printdiff('Cext',cext1,cext2)
    d += printdiff('Qext',qext1,qext2)
    d += printdiff('Cabs',cabs1,cabs2)
    d += printdiff('Qabs',qabs1,qabs2)
    fr.write('\tCrosSec-'+pol+': done')
    return d


def compare(line1,line2):
    
    adda_run(1,line1)
    adda_run(2,line2)
    
    print('\nCompare:\n'+line1+'\n'+line2)
    fr.write('\n\nCompare:\n'+line1+'\n'+line2)
    dtotal = 0
    dtotal += printsearch('X')
    dtotal += printsearch('Y')
    if (dtotal == 0):
        print('\033[32m', '\nPassed', '\033[0m', sep='') #green stdout
    else:
        print('\033[31m', '\nNot passed (see bb_results.txt)', '\033[0m', sep='') #red
    fr.write('\n\nDone\n___')
    
    shutil.rmtree('out1')
    shutil.rmtree('out2')

    
fname = "bb_results.txt"
fr = open(fname, "w")
fr.write("The comparison of equivalent cases (command lines) \nfor the scattering of Bessel beams calculated in ADDA.")
fr.write("\nDifferences are shown when |diff| > " + str(fdiff) + "\n\n___")
fr.close()

fr = open(fname, "a")


print('\n\nPlane-wave limit of LE Bessel beam')
fr.write('\n\nPlane-wave limit of LE Bessel beam')
opt1 = ' -sym no'
opt2 = ' -beam besselLE 0 0'
compare(opt1,opt2)

print('\n\nPlane-wave limit of LM Bessel beam')
fr.write('\n\nPlane-wave limit of LM Bessel beam')
opt1 = ' -sym no'
opt2 = ' -beam besselLM 0 0'
compare(opt1,opt2)

print('\n\nPlane-wave limit of CS Bessel beam')
fr.write('\n\nPlane-wave limit of CS Bessel beam')
opt1 = ' -sym no'
opt2 = ' -beam besselCS 0 0'
compare(opt1,opt2)

al1 = 5
al2 = 85

print('\n\nGeneralized and LE Bessel beams')
fr.write('\n\nGeneralized and LE Bessel beams')
opt1 = ' -beam besselM 2 '+str(al1)+' 0 0 0 1'
opt2 = ' -beam besselLE 2 '+str(al1)
compare(opt1,opt2)

opt1 = ' -beam besselM 2 '+str(al2)+' 0 0 0 1'
opt2 = ' -beam besselLE 2 '+str(al2)
compare(opt1,opt2)

print('\n\nGeneralized and LM Bessel beams')
fr.write('\n\nGeneralized and LM Bessel beams')
opt1 = ' -beam besselM 2 '+str(al1)+' 0 1 0 0'
opt2 = ' -beam besselLM 2 '+str(al1)
compare(opt1,opt2)

opt1 = ' -beam besselM 2 '+str(al2)+' 0 1 0 0'
opt2 = ' -beam besselLM 2 '+str(al2)
compare(opt1,opt2)

print('\n\nGeneralized and CS Bessel beams')
fr.write('\n\nGeneralized and CS Bessel beams')
opt1 = ' -beam besselM 2 '+str(al1)+' 0.5 0 0 0.5'
opt2 = ' -beam besselCS 2 '+str(al1)
compare(opt1,opt2)

opt1 = ' -beam besselM 2 '+str(al2)+' 0.5 0 0 0.5'
opt2 = ' -beam besselCS 2 '+str(al2)
compare(opt1,opt2)

print("\n\nGeneralized and CS' Bessel beams")
fr.write("\n\nGeneralized and CS' Bessel beams")
opt1 = ' -beam besselM 2 '+str(al1)+' 0.5 0 0 -0.5'
opt2 = ' -beam besselCSp 2 '+str(al1)
compare(opt1,opt2)

opt1 = ' -beam besselM 2 '+str(al2)+' 0.5 0 0 -0.5'
opt2 = ' -beam besselCSp 2 '+str(al2)
compare(opt1,opt2)

al1 = 9 #TEL an TML types have a bigger diff for smaller angles
al2 = 85

print('\n\nGeneralized and TEL Bessel beams')
fr.write('\n\nGeneralized and TEL Bessel beams')
opt1 = ' -beam besselM 2 '+str(al1)+' '+str(-1/math.sin(al1*math.pi/180))+' 0 0 '+str(1/math.tan(al1*math.pi/180))
#opt1 = ' -beam besselM 2 5 -11.4737132456698561 0 0 11.430052302761343'
opt2 = ' -beam besselTEL 2 '+str(al1)
compare(opt1,opt2)

opt1 = ' -beam besselM 2 '+str(al2)+' '+str(-1/math.sin(al2*math.pi/180))+' 0 0 '+str(1/math.tan(al2*math.pi/180))
opt2 = ' -beam besselTEL 2 '+str(al2)
compare(opt1,opt2)

print('\n\nGeneralized and TML Bessel beams')
fr.write('\n\nGeneralized and TML Bessel beams')
opt1 = ' -beam besselM 2 '+str(al1)+' 0 '+str(1/math.tan(al1*math.pi/180))+' '+str(1/math.sin(al1*math.pi/180))+' 0'
opt2 = ' -beam besselTML 2 '+str(al1)
compare(opt1,opt2)

opt1 = ' -beam besselM 2 '+str(al2)+' 0 '+str(1/math.tan(al2*math.pi/180))+' '+str(1/math.sin(al2*math.pi/180))+' 0'
opt2 = ' -beam besselTML 2 '+str(al2)
compare(opt1,opt2)

fr.close()
