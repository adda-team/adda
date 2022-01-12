from WKBr import find_wkb_ef
from rearrange_ef_file import rearrange_ef

# The program prepares the WKBr field.
# SET YOUR PARAMETERS:
size = 10
grid = 16
mi = 0
mlist = [1.01, 1.05, 1.1]
mytype = "complex_wkbr"
path = "./"

if __name__ == "__main__":
    for m in mlist:
        print("Current size = ", size)
        print("Current grid = ", grid)
        print("Current m = ", m)
        tail = str(size) + "-" + str(m) + "-" + str(grid) + ".dat"
        find_wkb_ef(0, 0, 0, m, mi, size / 2, 1, path + mytype + "-" + tail, grid, type=mytype, find_grid=True)

    # change columns to get y-polarization:
    for m in mlist:
        tail = str(size) + "-" + str(m) + "-" + str(grid)
        rearrange_ef(path + mytype + "-" + tail + ".dat", path + mytype + "-" + tail + "-Y-component.dat")