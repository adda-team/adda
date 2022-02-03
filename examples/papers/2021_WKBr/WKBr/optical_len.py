from iterative_method import iterative_method
from Region import region
from Fresnel_coefficients import transmission_coefficient, effective_refractive_indices
from Coordinate_systems import coordinates_in_meridional_plane, find_rotation_angle_in_mp, coordinates_in_lab_plane, \
    delta_angle, cos_refracted_angle_mp


# laboratory coordinate system - zero at the center of the sphere.
# x2_lab, y2_lab, z2_lab - coordinates of the point where the ray came in (lab. c.s.).
# z_ - the coordinate of the ray entrance into the sphere (WKB).
# type_ - "wkb+refraction", "analytic", "discrete".
# solution - "polynom", "approximate", "iterative".
def optical_len(x2_lab, y2_lab, z2_lab, z_, m, mi, type_, lines, k, solution="iterative"):
    # Go to the meridian plane:
    y2, z2 = coordinates_in_meridional_plane(x2_lab, y2_lab, z2_lab)
    # Find rotation angle:
    rotation_angle = find_rotation_angle_in_mp(x2_lab, y2_lab)
    # For WKBr:
    if type_ != "analytic" and type_ != "discrete":
        try:
            # Defined output parameters.
            info = "not_success"
            y1 = 0
            z1 = 0
            y1_2 = 0
            z1_2 = 0
            # We use an iterative algorithm:
            if solution == "iterative":
                if type_ == "wkb+refraction4" or type_ == "wkb+refraction5" or type_ == "wkb+refraction11":
                    # The first option: in the 1,2 solution region, we find only one root.
                    # In the 0 solution we use WKB.
                    cur_region = region(y2, z2, m, lines)
                    if cur_region == "one_root" or cur_region == "two_roots":
                        cur_region = "one_root"
                        info, y1, z1, y1_2, z1_2 = iterative_method(y2, z2, m, "one_root")
                    elif cur_region == "no_root":
                        info, y1, z1, y1_2, z1_2 = iterative_method(y2, z2, m, "no_root")
                    else:
                        print("Error in optical_len() function! In the first option.")
                elif type_ == "wkb+refraction" or type_ == "wkb+refraction6" or type_ == "wkb+refraction14" \
                        or type_ == "wkb+refraction12" or type_ == "wkb+refraction19"\
                        or type_ == "simple_wkbr":
                    # The first option: in the 1,2 solution region, we find only one root.
                    # In the 0 solution we use zero.
                    cur_region = region(y2, z2, m, lines)
                    if cur_region == "one_root" or cur_region == "two_roots":
                        cur_region = "one_root"
                        info, y1, z1, y1_2, z1_2 = iterative_method(y2, z2, m, "one_root")
                    elif cur_region == "no_root":
                        info = "no_calculation"
                    else:
                        print("Error in optical_len() function! In the first option.")
                elif type_ == "wkb+refraction15" or type_ == "wkb+refraction15-1" or type_ == "wkb+refraction15-2" \
                        or type_ == "wkb+refraction16" or type_ == "wkb+refraction16-1" or type_ == "wkb+refraction16-2" \
                        or type_ == "wkb+refraction17" or type_ == "wkb+refraction17-1" or type_ == "wkb+refraction17-2" \
                        or type_ == "wkb+refraction18" or type_ == "wkb+refraction18-1" or type_ == "wkb+refraction18-2" \
                        or type_ == "complex_wkbr":
                    # The 1-st option: R1 - use one solution. R2 - sum. R0 - use zero solution.
                    # The 3-rd option: R1 - one solution. R2 - rotation, transmission coefficients, sum. R0 - zero.
                    # The 4-rd option: R1 - one solution.
                    # R2 - rotation, transmission coefficients, convergence factor, sum. R0 - zero.
                    cur_region = region(y2, z2, m, lines)
                    if cur_region == "one_root" or cur_region == "two_roots":
                        info, y1, z1, y1_2, z1_2 = iterative_method(y2, z2, m, cur_region)
                    elif cur_region == "no_root":
                        info = "no_calculation"
                    else:
                        print("Error in optical_len() function! In wkb+refraction15(16,17).")
                elif type_ == "wkb+refraction2" or type_ == "wkb+refraction7" or type_ == "wkb+refraction13" \
                        or type_ == "wkb+refraction13-2":
                    # The second option: in the 1,2 solution region, we find one root and two roots respectively.
                    # In the 0 solution we use WKB.
                    cur_region = region(y2, z2, m, lines)
                    info, y1, z1, y1_2, z1_2 = iterative_method(y2, z2, m, cur_region)
                elif type_ == "wkb+refraction3" or type_ == "wkb+refraction10":
                    # The 3-rd option: in the 1,2 solution region, we find only one root.
                    # In the 0 solution we zero electric field.
                    cur_region = region(y2, z2, m, lines)
                    if cur_region == "one_root" or cur_region == "two_roots":
                        cur_region = "one_root"
                        info, y1, z1, y1_2, z1_2 = iterative_method(y2, z2, m, "one_root")
                    elif cur_region == "no_root":
                        info = "no_calculation"
                    else:
                        print("Error in optical_len() function! In the 3-rd option.")
                elif type_ == "wkb+refraction8" or type_ == "wkb+refraction9":
                    # The 4-th option: in the 1,2 solution region, we find one root and two roots respectively.
                    # In the 0 solution we zero electric field.
                    cur_region = region(y2, z2, m, lines)
                    if cur_region == "one_root" or cur_region == "two_roots":
                        info, y1, z1, y1_2, z1_2 = iterative_method(y2, z2, m, cur_region)
                    elif cur_region == "no_root":
                        info = "no_calculation"
                else:
                    print("Error in optical_len() function! No such type. ")
            else:
                print("Error in optical_len() function! In if solution == iterative")
        except RuntimeError as error:
            print("Oops! An error has occurred.")
            print(error)
            print("y2, z2 - ", y2, z2)
            info = "not_success"

        if info == "success":
            if cur_region == "one_root" or cur_region == "no_root":
                l1 = abs(z1 + 1)
                l2 = ((y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5
                l1_2 = 0
                l2_2 = 0
                # Nr, Ni = find_adjusted_refractive_indices(m, mi, y1 / radius)
                t_per1, t_par1 = transmission_coefficient(y1, m, mi, k)
                t_per2, t_par2 = 0, 0
                x1_lab, y1_lab, z1_lab = coordinates_in_lab_plane(y1, z1, rotation_angle)
                da1 = delta_angle(x1_lab, m)
                da2 = 0
                cos_t1 = cos_refracted_angle_mp(y1, m)
                cos_t2 = 0
                N1, K1 = effective_refractive_indices(m, mi, cos_t1)
                N2, K2 = 0, 0
            elif cur_region == "two_roots":
                l1 = abs(z1 + 1)
                l2 = ((y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5
                l1_2 = abs(z1_2 + 1)
                l2_2 = ((y2 - y1_2) ** 2 + (z2 - z1_2) ** 2) ** 0.5
                # Nr, Ni = find_adjusted_refractive_indices(m, mi, y1 / radius)
                # Nr2, Ni2 = find_adjusted_refractive_indices(m, mi, y1_2 / radius)
                t_per1, t_par1 = transmission_coefficient(y1, m, mi, k)
                t_per2, t_par2 = transmission_coefficient(y1_2, m, mi, k)
                x1_lab, y1_lab, z1_lab = coordinates_in_lab_plane(y1, z1, rotation_angle)
                x1_lab_2, y1_lab_2, z1_lab_2 = coordinates_in_lab_plane(y1_2, z1_2, rotation_angle)
                da1 = delta_angle(x1_lab, m)
                da2 = delta_angle(x1_lab_2, m)
                cos_t1 = cos_refracted_angle_mp(y1, m)
                cos_t2 = cos_refracted_angle_mp(y1_2, m)
                N1, K1 = effective_refractive_indices(m, mi, cos_t1)
                N2, K2 = effective_refractive_indices(m, mi, cos_t2)
            else:
                print("Error in optical_len() function!")
        elif info == "not_success":
            type_ = "analytic"
        elif info == "no_calculation":
            # This code is executed if we do not want to calculate anything and zero the electric field.
            l1, l2, l1_2, l2_2 = 0, 0, 0, 0
            t_per1, t_per2, t_par1, t_par2 = 0, 0, 0, 0
            da1, da2, cos_t1, cos_t2 = 0, 0, 0, 0
            N1, K1, N2, K2 = 0, 0, 0, 0
        else:
            print("Error in optical_len() function! Unknown parameter info. ")
    # For WKB:
    if type_ == "analytic" or type_ == "discrete":
        l1 = abs(z_ + 1)
        l2 = abs(z2_lab - z_)
        l1_2 = 0
        l2_2 = 0
        cur_region = "one_root"
        # In the case of WKB, there is no refraction, so instead of y1 use y2.
        t_per1, t_par1 = transmission_coefficient(y2, m, mi, k)
        t_per2, t_par2 = t_per1, t_par1
        x1_lab, y1_lab, z1_lab = coordinates_in_lab_plane(y2, z2, rotation_angle)
        da1 = delta_angle(x1_lab, m)
        da2 = da1
        cos_t1 = cos_refracted_angle_mp(y2, m)
        cos_t2 = cos_refracted_angle_mp(y2, m)
        N1, K1 = effective_refractive_indices(m, mi, cos_t1)
        N2, K2 = effective_refractive_indices(m, mi, cos_t2)

    return l1, l2, l1_2, l2_2, type_, cur_region, t_per1, t_per2, t_par1, t_par2, da1, da2, rotation_angle, \
           cos_t1, cos_t2, N1, K1, N2, K2
