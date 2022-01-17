from sys import argv
from WKBr import find_wkb_ef
from rearrange_ef_file import rearrange_ef

# The program prepares the WKBr field.
# launch example: python main.py 10 1.1 16 simple_wkbr

if __name__ == "__main__":
    size = float(argv[1])
    m = float(argv[2])
    grid = float(argv[3])
    mytype = argv[4] # supported options: "simple_wkbr" (WKBrI), "complex_wkbr" (WKBrII), and "analytic" (conventional WKB)  
    print("Current size = ", size)
    print("Current grid = ", grid)
    print("Current m = ", m)
    print("Current WKBr type = ", mytype)
    tail = str(size) + "-" + str(m) + "-" + str(grid) + ".dat"
    find_wkb_ef(0, 0, 0, m, 0, size / 2, 1, mytype + "-" + tail, grid, type=mytype, find_grid=True)

    # change columns to get y-polarization:
    tail = str(size) + "-" + str(m) + "-" + str(grid)
    rearrange_ef(mytype + "-" + tail + ".dat", mytype + "-" + tail + "-Y-component.dat")