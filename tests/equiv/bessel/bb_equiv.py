'''
# This code compares equivalent cases (command lines)
# for the scattering of Bessel beams calculated in ADDA.
'''

import os,shutil,re,math


# path to adda executable
#adda_exec = "../../win64/adda.exe"
adda_exec = os.path.abspath(__file__ + "/../../../../src/seq/adda")

dirname = 'out'
fdiff = 1e-6 # fixed difference

pi = 3.14159265359
                                                                              

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
    cmdline = adda_exec + ' -dir ' + dname + option
    os.system(cmdline)


# print difference
def printdiff(val,x,y):
    cdiff = math.fabs(x-y) # calculated difference
    if (cdiff > fdiff):
        print('\n\t\t'+val + ':\n\t\tcase 1:\t'+str(x)+'\n\t\tcase 2:\t'+str(y)+'\n\t\tdiff:\t'+str(cdiff)+'\n')
        #print('\ncase 1:\t'+str(x)+'\ncase 2:\t'+str(y)+'\nvariation:\t'+str(cdiff)+' %')
        fr.write('\n\t\t'+val + ':\n\t\tcase 1:\t'+str(x)+'\n\t\tcase 2:\t'+str(y)+'\n\t\tdiff:\t'+str(cdiff)+'\n')


def printsearch(pol):
    print('\n\tCrosSec-'+pol+': search diff...')
    fr.write('\n\n\tCrosSec-'+pol+': search diff')
    cext1,qext1,cabs1,qabs1 = extractCrosSeq(1,pol)
    cext2,qext2,cabs2,qabs2 = extractCrosSeq(2,pol)
    printdiff('Cext',cext1,cext2)
    printdiff('Qext',qext1,qext2)
    printdiff('Cabs',cabs1,cabs2)
    printdiff('Qabs',qabs1,qabs2)
    print('\tCrosSec-'+pol+': done')
    fr.write('\tCrosSec-'+pol+': done')


def compare(line1,line2):
    
    adda_run(1,line1)
    adda_run(2,line2)
    
    print('\nCompare:\n'+line1+'\n'+line2)
    fr.write('\n\nCompare:\n'+line1+'\n'+line2)
    printsearch('X')
    printsearch('Y')
    print('\nDone\n___')
    fr.write('\n\nDone\n___')
    
    shutil.rmtree('out1')
    shutil.rmtree('out2')

    
fr = open("results.txt", "w")
fr.write("The comparison of equivalent cases (command lines) \nfor the scattering of Bessel beams calculated in ADDA.")
fr.write("\nDifferences are shown when |diff| > " + str(fdiff) + "\n\n___")
fr.close()

fr = open("results.txt", "a")


print('\n\nPlane-wave limit of LE Bessel beam')
fr.write('\n\nPlane-wave limit of LE Bessel beam')
common = ' -sym no'
opt1 = common
opt2 = common + ' -beam besselLE 0 0'
compare(opt1,opt2)

print('\n\nPlane-wave limit of LM Bessel beam')
fr.write('\n\nPlane-wave limit of LM Bessel beam')
common = ' -sym no'
opt1 = common
opt2 = common + ' -beam besselLM 0 0'
compare(opt1,opt2)

print('\n\nPlane-wave limit of CS Bessel beam')
fr.write('\n\nPlane-wave limit of CS Bessel beam')
common = ' -sym no'
opt1 = common
opt2 = common + ' -beam besselCS 0 0'
compare(opt1,opt2)

print('\n\nGeneralized and LE Bessel beams')
fr.write('\n\nGeneralized and LE Bessel beams')
opt1 = ' -beam besselM 2 5 0 0 0 1'
opt2 = ' -beam besselLE 2 5'
compare(opt1,opt2)

opt1 = ' -beam besselM 2 85 0 0 0 1'
opt2 = ' -beam besselLE 2 85'
compare(opt1,opt2)

print('\n\nGeneralized and LM Bessel beams')
fr.write('\n\nGeneralized and LM Bessel beams')
opt1 = ' -beam besselM 2 5 0 1 0 0'
opt2 = ' -beam besselLM 2 5'
compare(opt1,opt2)

opt1 = ' -beam besselM 2 85 0 1 0 0'
opt2 = ' -beam besselLM 2 85'
compare(opt1,opt2)

print('\n\nGeneralized and CS Bessel beams')
fr.write('\n\nGeneralized and CS Bessel beams')
opt1 = ' -beam besselM 2 5 0.5 0 0 0.5'
opt2 = ' -beam besselCS 2 5'
compare(opt1,opt2)

opt1 = ' -beam besselM 2 85 0.5 0 0 0.5'
opt2 = ' -beam besselCS 2 85'
compare(opt1,opt2)

print("\n\nGeneralized and CS' Bessel beams")
fr.write("\n\nGeneralized and CS' Bessel beams")
opt1 = ' -beam besselM 2 5 0.5 0 0 -0.5'
opt2 = ' -beam besselCSp 2 5'
compare(opt1,opt2)

opt1 = ' -beam besselM 2 85 0.5 0 0 -0.5'
opt2 = ' -beam besselCSp 2 85'
compare(opt1,opt2)

print('\n\nGeneralized and TEL Bessel beams')
fr.write('\n\nGeneralized and TEL Bessel beams')
opt1 = ' -beam besselM 2 5 '+str(-1/math.sin(5*pi/180))+' 0 0 '+str(1/math.tan(5*pi/180))
opt2 = ' -beam besselTEL 2 5'
compare(opt1,opt2)

opt1 = ' -beam besselM 2 85 '+str(-1/math.sin(85*pi/180))+' 0 0 '+str(1/math.tan(85*pi/180))
opt2 = ' -beam besselTEL 2 85'
compare(opt1,opt2)

print('\n\nGeneralized and TML Bessel beams')
fr.write('\n\nGeneralized and TML Bessel beams')
opt1 = ' -beam besselM 2 5 0 '+str(1/math.tan(5*pi/180))+' '+str(1/math.sin(5*pi/180))+' 0'
opt2 = ' -beam besselTML 2 5'
compare(opt1,opt2)

opt1 = ' -beam besselM 2 85 0 '+str(1/math.tan(85*pi/180))+' '+str(1/math.sin(85*pi/180))+' 0'
opt2 = ' -beam besselTML 2 85'
compare(opt1,opt2)

fr.close()