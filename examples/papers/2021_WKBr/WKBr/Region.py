from math import cos, asin, tan
from iterative_method import func, func2


# This class stores the coefficients of the lines that bound regions with one solution, with two and without solutions.
# Radius = 1.
# z = a * y + b
class BoundaryLines(object):
    def __init__(self, m):
        self.m = m
        # grazing ray:
        tmp = 1 / (m ** 2 - 1) ** 0.5
        self.a1 = - tmp
        self.b1 = tmp
        # Descartes ray:
        sin_tmp = ((4 - m ** 2) / 3) ** 0.5
        teta = asin(sin_tmp)
        psi = 2 * asin(sin_tmp / m) - teta
        k1 = tan((teta - psi) / 2)
        self.a2 = - 1 / k1
        self.b2 = sin_tmp / k1 - cos(teta)


# check that the point in the sphere.
# R = 1.
def in_sphere(y, z):
    r2 = y ** 2 + z ** 2
    if r2 <= 1:
        return True
    else:
        return False


# Function for determining in which region a given point is located.
# R = 1.
# Returns: "one_root", "two_roots", "no_root", "error"
def region(y2, z2, m, lines):
    # empirical Descartes ray thickness parameter:
    delta = 0.01 + (m - 1) / 10  # 0.05
    if in_sphere(y2, z2):
        if z2 >= 0:
            if m == lines.m:
                z_l1 = lines.a1 * y2 + lines.b1  # grazing ray
                y2_l2 = (z2 - lines.b2) / lines.a2  # Descartes ray
            else:
                print("Error in region() function!")
                return "error"
            if z2 >= z_l1:
                if y2 < y2_l2 - delta:
                    return "two_roots"
                elif y2 > y2_l2 + delta:
                    return "no_root"
                elif y2_l2 - delta <= y2 <= y2_l2 + delta:
                    return region2(y2, z2, m)
            else:
                return "one_root"
        else:
            return "one_root"
    else:
        print("Error in region() function! Point is not in sphere. Current (y2, z2) = ", y2, z2)
        return "error"


def in_square(y, z, R):
    if 0 <= y <= R and -R <= z <= 0:
        return True
    else:
        return False


# Function for determining in which region a given point is located (without grazing and Descartes rays).
# Returns: "one_root", "two_roots", "no_root", "error".
def region2(y2, z2, m):
    if in_sphere(y2, z2):
        if z2 > 0:
            root_flag, y1, z1 = func(y2, z2, m, y2)
            if root_flag:
                # f(1) > g(1)
                f_1 = (1 - y2) / z2
                g_1 = (m ** 2 - 1) ** 0.5
                if f_1 > g_1:
                    return "one_root"
                else:
                    root_flag2, y1_2, z1_2 = func2(y2, z2, m, 1)
                    if root_flag2:
                        return "two_roots"
                    else:
                        return "one_root"
            else:
                return "no_root"
        else:
            return "one_root"
    else:
        print("Error in region() function! Point is not in sphere. Current (y2, z2) = ", y2, z2)
        return "error"

