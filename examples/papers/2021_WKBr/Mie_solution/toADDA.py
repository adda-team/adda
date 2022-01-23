# This script casts grid Bhfield / Scattnlay to the ADDA type
# launch example: python toADDA.py -s 10 -m 1.1 -g 16 -t scattnlay

import sys, getopt
import linecache
from math import isnan

# path1 - path to Scattnlay or Bhfield file.
# path2 - path to the file where to save.
# type - "bhfield", "scattnlay".
# R - sphere radius.
def scattnlay_bhfield_to_adda(path1, path2, type, R, grid):
    count_nan = 0

    f1 = open(path1)
    if type == "bhfield":
        f1.readline()
        f1.readline()
    elif type != "scattnlay" and type != "bhfield":
        print("Error in scattnlay_bhfield_to_adda()! Unknown type.")

    f2 = open(path2, 'w')
    f2.write('x y z |E|^2 Ex.r Ex.i Ey.r Ey.i Ez.r Ez.i \n')

    if type == "bhfield":
        while True:
            numbers1 = f1.readline()
            if len(numbers1) == 0:  # EOF
                break
            numbers1 = numbers1.split()
            if len(numbers1) <= 1:
                break
            x1 = float(numbers1[0])
            y1 = float(numbers1[1])
            z1 = float(numbers1[2])
            R1 = (x1 ** 2 + y1 ** 2 + z1 ** 2) ** 0.5
            if R1 > R:
                continue

            exr = float(numbers1[3])
            eyr = float(numbers1[4])
            ezr = float(numbers1[5])
            exi = float(numbers1[6])
            eyi = float(numbers1[7])
            ezi = float(numbers1[8])
            e = exr ** 2 + exi ** 2 + eyr ** 2 + eyi ** 2 + ezr ** 2 + ezi ** 2
            f2.write(str(x1) + ' ' + str(y1) + ' ' + str(z1) + ' ' +
                     str(e) + ' ' +
                     str(exr) + ' ' + str(exi) + ' ' +
                     str(eyr) + ' ' + str(eyi) + ' ' +
                     str(ezr) + ' ' + str(ezi) + ' ' +
                     '\n')
    elif type == "scattnlay":
        for k in range(grid):
            for j in range(grid):
                for i in range(grid):
                    # 6 - number of lines to skip. 
                    # It is this number of lines (without data) that is observed in the Scattnlay version dated June 7, 2021.
                    numbers1 = linecache.getline(path1, 6 + 1 + grid * grid * i + grid * j + k)

                    if len(numbers1) == 0:  # EOF
                        break

                    if len(numbers1) <= 1: # To remove an extra line in scattnlay at the end of the file.
                        break
                    numbers1 = numbers1.split(", ")

                    if len(numbers1) <= 1:
                        break

                    x1 = float(numbers1[0])
                    y1 = float(numbers1[1])
                    z1 = float(numbers1[2])
                    R1 = (x1 ** 2 + y1 ** 2 + z1 ** 2) ** 0.5
                    if R1 > R:
                        continue

                    exr = float(numbers1[3])
                    exi = float(numbers1[4])
                    eyr = float(numbers1[5])
                    eyi = float(numbers1[6])
                    ezr = float(numbers1[7])
                    ezi = float(numbers1[8])

                    if isnan(exr) or isnan(eyr) or isnan(ezr) or isnan(exi) or isnan(eyi) or isnan(ezi):
                        count_nan += 1
                        continue

                    e = exr ** 2 + exi ** 2 + eyr ** 2 + eyi ** 2 + ezr ** 2 + ezi ** 2
                    f2.write(str(x1) + ' ' + str(y1) + ' ' + str(z1) + ' ' +
                             str(e) + ' ' +
                             str(exr) + ' ' + str(exi) + ' ' +
                             str(eyr) + ' ' + str(eyi) + ' ' +
                             str(ezr) + ' ' + str(ezi) + ' ' +
                             '\n')
    else:
        print("Error in scattnlay_bhfield_to_adda()! Unknown type.")

    if count_nan != 0:
        print("There is nan! Number of nan: ", count_nan)
    f1.close()
    f2.close()


def main(argv):
    # default options:
    size = "10"
    m = "1.1"
    grid = "16"
    type = "scattnlay" # supported options: "scattnlay", "bhfield"
    # read options from command line:
    try:
        opts, args = getopt.getopt(argv,"hs:m:g:t:")
    except getopt.GetoptError:
        print ('python toADDA.py -s <size> -m <m> -g <grid> -t <type>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('python toADDA.py -s <size> -m <m> -g <grid> -t <type>')
            sys.exit()
        elif opt == "-s":
         size = arg
        elif opt == "-m":
         m = arg
        elif opt == "-g":
         grid = arg
        elif opt == "-t":
         type = arg 
    print("Current size = ", size)
    print("Current grid = ", grid)
    print("Current m = ", m)
    print("Current type = ", type)
    path_exact = type + "-" + size + "-" + m + "-" + grid + ".dat"
    path_exact_adda_x_component = "if-" + size + "-" + m + "-" + grid + "-X-component.dat"
    scattnlay_bhfield_to_adda(path_exact, path_exact_adda_x_component, type, float(size) / 2, int(grid))
    
    
if __name__ == '__main__':
    main(sys.argv[1:])

