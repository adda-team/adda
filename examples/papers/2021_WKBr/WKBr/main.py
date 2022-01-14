from WKBr import find_wkb_ef
from rearrange_ef_file import rearrange_ef

# The program prepares the WKBr field.
# SET YOUR PARAMETERS:
# TODO: should accept parameters through the command line (not sure about mi)
size = 10
m = 1.1
mi = 0
grid = 16

# Three options are supported: "simple_wkbr" (WKBrI), "complex_wkbr" (WKBrII), and "analytic" (conventional WKB)  
mytype = "complex_wkbr"
path = "./"

if __name__ == "__main__":
    print("Current size = ", size)
    print("Current grid = ", grid)
    print("Current m = ", m)
    tail = str(size) + "-" + str(m) + "-" + str(grid) + ".dat"
    find_wkb_ef(0, 0, 0, m, mi, size / 2, 1, path + mytype + "-" + tail, grid, type=mytype, find_grid=True)

    # change columns to get y-polarization:
    tail = str(size) + "-" + str(m) + "-" + str(grid)
    rearrange_ef(path + mytype + "-" + tail + ".dat", path + mytype + "-" + tail + "-Y-component.dat")