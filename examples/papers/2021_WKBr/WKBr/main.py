import sys, getopt
from WKBr import find_wkb_ef
from rearrange_ef_file import rearrange_ef

# The program prepares the WKBr field.
# launch example: python main.py -s 10 -m 1.1 -g 16 -t simple_wkbr

def main(argv):
    # default options:
    size = 10
    m = 1.1
    grid = 16
    type = "simple_wkbr" # supported options: "simple_wkbr" (WKBrI), "complex_wkbr" (WKBrII), and "analytic" (conventional WKB)
    # read options from command line:
    try:
        opts, args = getopt.getopt(argv,"hs:m:g:t:")
    except getopt.GetoptError:
        print ('python main.py -s <size> -m <m> -g <grid> -t <type>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('python main.py -s <size> -m <m> -g <grid> -t <type>')
            sys.exit()
        elif opt == "-s":
         size = float(arg)
        elif opt == "-m":
         m = float(arg)
        elif opt == "-g":
         grid = int(arg)
        elif opt == "-t":
         type = arg 
         
    print("Current size = ", size)
    print("Current grid = ", grid)
    print("Current m = ", m)
    print("Current WKBr type = ", type)
    tail = str(size) + "-" + str(m) + "-" + str(grid) + ".dat"
    find_wkb_ef(0, 0, 0, m, 0, size / 2, 1, type + "-" + tail, grid, type=type, find_grid=True)

    # change columns to get y-polarization:
    tail = str(size) + "-" + str(m) + "-" + str(grid)
    rearrange_ef(type + "-" + tail + ".dat", type + "-" + tail + "-Y-component.dat")
    
if __name__ == '__main__':
    main(sys.argv[1:])