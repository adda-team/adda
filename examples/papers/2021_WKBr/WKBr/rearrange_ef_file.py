from math import isnan


# Function for permutation of electric field components.
# Тип: ADDA.
# EX and EY are swapped.
def rearrange_ef(open_path, save_path):
    print("rearrange function is working")
    f1 = open(open_path, 'r')
    f2 = open(save_path, 'w')
    f2.write(f1.readline())
    while True:
        numbers = f1.readline()
        if len(numbers) == 0:  # EOF
            break
        numbers = numbers.split()
        x = float(numbers[0])
        y = float(numbers[1])
        z = float(numbers[2])
        e = float(numbers[3])
        exr = float(numbers[4])
        exi = float(numbers[5])
        s = str(x)+' '+str(y)+' '+str(z)+' '+str(e)+' 0 0 '+str(exr)+' '+str(exi)+' 0 0 \n'
        f2.write(s)
    f1.close()
    f2.close()


# The function corrects file f and saves to f2.
def correct_scattnlay_file(path, path2):
    f = open(path)
    f2 = open(path2, 'w')
    f2.write('x y z |E|^2 Ex.r Ex.i Ey.r Ey.i Ez.r Ez.i \n')
    counter = 0
    counter_nan = 0
    while True:
        counter += 1
        numbers_line = f.readline()
        if numbers_line[0:7] == "Warning" or numbers_line[0:11] == "         X,":
            continue
        numbers = numbers_line.split(", ")
        if len(numbers) <= 1:
            break
        x = float(numbers[0])
        y = float(numbers[1])
        z = float(numbers[2])
        exr = float(numbers[3])
        exi = float(numbers[4])
        eyr = float(numbers[5])
        eyi = float(numbers[6])
        ezr = float(numbers[7])
        ezi = float(numbers[8])
        if isnan(exr) or isnan(eyr) or isnan(ezr) or isnan(exi) or isnan(eyi) or isnan(ezi):
            counter_nan += 1
            exr = 0
            exi = 0
            eyr = 0
            eyi = 0
            ezr = 0
            ezi = 0
            f2.write(str(x) + ', ' + str(y) + ', ' + str(z) + ', ' +
                     str(exr) + ', ' + str(exi) + ', ' + str(eyr) + ', ' +
                     str(eyi) + ', ' + str(ezr) + ', ' + str(ezi) + ' ' +
                     "\n")
        else:
            f2.write(numbers_line)
    f.close()
    f2.close()
    print("Total number in cube: ", counter)
    print("Total number in sphere: ", 3*counter/(3.1415*2**0.5))
    print("nan number: ", counter_nan)
