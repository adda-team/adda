import sys

if __name__ == "__main__":
    size = float(sys.argv[1])
    grid = float(sys.argv[2])
    
    dip = size / grid
    b = size / 2. - dip / 2.
    
    lb = str(-b)
    rb = str(b)
    print(rb)


