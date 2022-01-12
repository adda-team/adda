import numpy as np
from math import sin, cos, exp, acos, asin, pi
from Region import BoundaryLines
from optical_len import optical_len
from Coordinate_systems import apply_rotation_electric_field_vector, sum_ef, apply_transmission_coefficient


# Function for finding the argument of the exponent of the electric field of the transmitted wave:
# l1 - outside the particle
# l2 - in the particle
def find_arg(k, radius, m, l1, l2):
    return k * radius * (-1 + l1 + m * l2)


# Function for finding the attenuation of the electric field of the transmitted wave:
def find_attenuation(k, l2, mi, cos_t, radius):
    return exp(- k * mi * l2 * radius * cos_t)


def safe_div(x, y):
    if abs(y) < 1e-16:
        print("Error: division by zero")
        return 0
    return x / y


# Function for finding the "convergence factor" due to the curvature.
# See details: Graeme L. James "Geometrical Theory of Diffraction for Electromagnetic Waves" (1986).
# l - optical path in the particle
def find_convergence_factor(radius, m, l, inc_angle, refracted_angle):
    C = np.array([[1 / radius, 0], [0, 1 / radius]])  # surface curvature matrix
    # Qi = np.array([[0, 0], [0, 0]])  # incident front curvature matrix
    # Ki = np.array([[1, 0], [0, -cos(inc_angle)]])  # matrix transverse coordinates
    RevKt = np.array([[1, 0], [0, -1 / cos(refracted_angle)]])
    Gi = C * cos(inc_angle)  # + Ki*Qi*Ki
    Gt = Gi / m
    Qt = RevKt.dot(Gt - C * cos(refracted_angle)).dot(RevKt)
    Q11 = Qt[0][0]
    Q22 = Qt[1][1]
    Q12 = Qt[0][1]
    den1 = Q11 + Q22 + ((Q11 - Q22) ** 2 + 4 * Q12 ** 2) ** 0.5
    den2 = Q11 + Q22 - ((Q11 - Q22) ** 2 + 4 * Q12 ** 2) ** 0.5
    ro1 = safe_div(2, den1)
    ro2 = safe_div(2, den2)
    K = safe_div(ro1 * ro2, (ro1 + l) * (ro2 + l))
    if K < 0:
        negative_sign = True
    else:
        negative_sign = False
    K = abs(K)
    K = K ** 0.5
    return K, negative_sign


# Function for applying the convergence factor.
def apply_convergence_factor(exr_new, exi_new, eyr_new, eyi_new, ezr_new, ezi_new, K):
    exr_new *= K
    exi_new *= K
    eyr_new *= K
    eyi_new *= K
    ezr_new *= K
    ezi_new *= K
    return exr_new, exi_new, eyr_new, eyi_new, ezr_new, ezi_new


# A plane wave is incident along the z axis.
# The wave is linearly polarized along the x axis.
def find_wkb_ef(x_arr, y_arr, z_arr, m, mi, radius, k, path, grid, type="analytic", find_grid=False,
                solution_method="iterative"):
    if type == "analytic" or type == "discrete":
        print("WKB calculation: " + type)
    else:
        print("WKBr calculation: " + type)
    # Find electric field WKB:
    f = open(path, 'w')
    f.write('x y z |E| Ex.r Ex.i Ey.r Ey.i Ez.r Ez.i \n')
    num_one_root = 0
    num_two_roots = 0
    num_no_roots = 0

    diameter = 2 * radius
    d = diameter / grid
    start = (-diameter + d) / 2
    finish = (diameter - d) / 2
    i = 0
    x_cur = start - d
    y_cur = start
    z_cur = start
    while True:
        if find_grid:
            # Build a coordinate grid.
            x_cur += d
            if x_cur > finish:
                x_cur = start
                y_cur += d
                if y_cur > finish:
                    print('z = ', z_cur)
                    y_cur = start
                    z_cur += d
                    if z_cur > finish:
                        break
            R_cur = (x_cur ** 2 + y_cur ** 2 + z_cur ** 2) ** 0.5
            if radius < R_cur:
                continue
        else:
            # Use the already made grid.
            if i == len(x_arr):
                # End of calculations:
                break
            else:
                x_cur = x_arr[i]
                y_cur = y_arr[i]
                z_cur = z_arr[i]
                i += 1

        # Normalize coordinates to get dimensionless:
        x_cur_old = x_cur
        y_cur_old = y_cur
        z_cur_old = z_cur
        x_cur /= radius
        y_cur /= radius
        z_cur /= radius
        # Separate regions using rays:
        lines = BoundaryLines(m)

        # The boundary is set analytically (not like in ADDA).
        # The usual sphere equation is used here. 0 coordinate system at the center of the sphere.
        z_ = - (1 - x_cur ** 2 - y_cur ** 2) ** 0.5
        l1, l2, l1_2, l2_2, type2, cur_region, t_per1, t_per2, t_par1, t_par2, da1, da2, rotation_angle, \
        cos_t1, cos_t2, N1, K1, N2, K2 = optical_len(x_cur,
                                                     y_cur,
                                                     z_cur,
                                                     z_,
                                                     m,
                                                     mi,
                                                     type,
                                                     lines,
                                                     k,
                                                     solution=solution_method)

        if type == "analytic" or type == "discrete":
            arg = find_arg(k, radius, N1, l1, l2)
            attenuation1 = find_attenuation(k, l2, K1, cos_t1, radius)
            exr_new = cos(arg) * attenuation1
            exi_new = sin(arg) * attenuation1
            eyr_new, eyi_new, ezr_new, ezi_new = 0, 0, 0, 0
            e = (exr_new ** 2 + exi_new ** 2 + eyr_new ** 2 + eyi_new ** 2 + ezr_new ** 2 + ezi_new ** 2) ** 0.5
            num_one_root += 1
        elif type == "simple_wkbr":
            # WKBr with rotation of the electric field.
            if cur_region == "one_root" or cur_region == "two_roots":
                arg = find_arg(k, radius, N1, l1, l2)
                attenuation1 = find_attenuation(k, l2, K1, cos_t1, radius)
                exr_new = cos(arg) * attenuation1
                exi_new = sin(arg) * attenuation1
                exr_new, exi_new, eyr_new, eyi_new, ezr_new, ezi_new = \
                    apply_rotation_electric_field_vector(exr_new, exi_new, da1)
                e = (exr_new ** 2 + exi_new ** 2 + eyr_new ** 2 + eyi_new ** 2 + ezr_new ** 2 + ezi_new ** 2) ** 0.5
                num_one_root += 1
            elif cur_region == "no_root":
                exr_new, exi_new, eyr_new, eyi_new, ezr_new, ezi_new, e = 0, 0, 0, 0, 0, 0, 0
                num_no_roots += 1
        elif type == "complex_wkbr":
            if cur_region == "one_root":
                arg = find_arg(k, radius, N1, l1, l2)
                attenuation1 = find_attenuation(k, l2, K1, cos_t1, radius)
                exr_new = cos(arg) * attenuation1
                exi_new = sin(arg) * attenuation1
                exr_new, exi_new, eyr_new, eyi_new = \
                    apply_transmission_coefficient(exr_new, exi_new, t_per1, t_par1, rotation_angle)
                exr_new, exi_new, eyr_new, eyi_new, ezr_new, ezi_new = \
                    apply_rotation_electric_field_vector(exr_new, exi_new, da1)
                ref_ang = acos(cos_t1)
                inc_ang = asin(m * sin(ref_ang))
                K, negative_sign = find_convergence_factor(radius, m, radius * l2, inc_ang, ref_ang)
                if negative_sign:
                    print("Negative sign in one-solution region!")
                exr_new, exi_new, eyr_new, eyi_new, ezr_new, ezi_new = \
                    apply_convergence_factor(exr_new, exi_new, eyr_new, eyi_new, ezr_new, ezi_new, K)
                e = (exr_new ** 2 + exi_new ** 2 + eyr_new ** 2 + eyi_new ** 2 + ezr_new ** 2 + ezi_new ** 2) ** 0.5
                num_one_root += 1
            elif cur_region == "two_roots":
                attenuation1 = find_attenuation(k, l2, K1, cos_t1, radius)
                attenuation2 = find_attenuation(k, l2_2, K2, cos_t2, radius)
                arg1 = find_arg(k, radius, N1, l1, l2)
                arg2 = find_arg(k, radius, N2, l1_2, l2_2)

                ref_ang1 = acos(cos_t1)
                inc_ang1 = asin(m * sin(ref_ang1))
                K1, negative_sign1 = find_convergence_factor(radius, m, radius * l2, inc_ang1, ref_ang1)
                if negative_sign1:
                    print("Negative sign in double-solution region (first ray)!")

                ref_ang2 = acos(cos_t2)
                inc_ang2 = asin(m * sin(ref_ang2))
                K2, negative_sign2 = find_convergence_factor(radius, m, radius * l2_2, inc_ang2, ref_ang2)
                if negative_sign2:
                    arg2 -= pi / 2

                exr_new1 = cos(arg1) * attenuation1
                exr_new2 = cos(arg2) * attenuation2
                exi_new1 = sin(arg1) * attenuation1
                exi_new2 = sin(arg2) * attenuation2
                # -----------------------------------------------------------#
                # First ray
                exr_new1, exi_new1, eyr_new1, eyi_new1 = \
                    apply_transmission_coefficient(exr_new1, exi_new1, t_per1, t_par1, rotation_angle)
                exr_new1, exi_new1, eyr_new1, eyi_new1, ezr_new1, ezi_new1 = \
                    apply_rotation_electric_field_vector(exr_new1, exi_new1, da1)
                exr_new1, exi_new1, eyr_new1, eyi_new1, ezr_new1, ezi_new1 = \
                    apply_convergence_factor(exr_new1, exi_new1, eyr_new1, eyi_new1, ezr_new1, ezi_new1, K1)
                # -----------------------------------------------------------#
                # Second ray
                exr_new2, exi_new2, eyr_new2, eyi_new2 = \
                    apply_transmission_coefficient(exr_new2, exi_new2, t_per2, t_par2, rotation_angle)
                exr_new2, exi_new2, eyr_new2, eyi_new2, ezr_new2, ezi_new2 = \
                    apply_rotation_electric_field_vector(exr_new2, exi_new2, da2)
                exr_new2, exi_new2, eyr_new2, eyi_new2, ezr_new2, ezi_new2 = \
                    apply_convergence_factor(exr_new2, exi_new2, eyr_new2, eyi_new2, ezr_new2, ezi_new2, K2)
                # -----------------------------------------------------------#
                exr_new, exi_new, eyr_new, eyi_new, ezr_new, ezi_new = sum_ef(exr_new1, exr_new2, exi_new1, exi_new2,
                                                                              eyr_new1, eyr_new2, eyi_new1, eyi_new2,
                                                                              ezr_new1, ezr_new2, ezi_new1, ezi_new2)
                e = (exr_new ** 2 + exi_new ** 2 + eyr_new ** 2 + eyi_new ** 2 + ezr_new ** 2 + ezi_new ** 2) ** 0.5
                num_two_roots += 1
            elif cur_region == "no_root":
                exr_new, exi_new, eyr_new, eyi_new, ezr_new, ezi_new, e = 0, 0, 0, 0, 0, 0, 0
                num_no_roots += 1
        else:
            print("Error in find_wkb_ef() function!")
        # Write electric field
        f.write(str(x_cur_old) + ' ' + str(y_cur_old) + ' ' + str(z_cur_old) + ' ' + str(e) + ' '
                + str(exr_new) + ' ' + str(exi_new) + ' '
                + str(eyr_new) + ' ' + str(eyi_new) + ' '
                + str(ezr_new) + ' ' + str(ezi_new) + '\n')
        # you need to return the values, otherwise an error occurs when generating the grid:
        x_cur = x_cur_old
        y_cur = y_cur_old
        z_cur = z_cur_old
    f.close()

    print("Number of dots in one root region = ", num_one_root)
    print("Number of dots in two root region = ", num_two_roots)
    print("Number of dots in no root region = ", num_no_roots)
    sum_ = num_no_roots + num_one_root + num_two_roots
    print("Number of all dots = ", sum_)

    return sum_, num_no_roots, num_one_root, num_two_roots
