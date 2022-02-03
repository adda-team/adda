tolerance = 1e-8
max_iterations = 1500


# Exact solution using iterative method.
# The function returns the coordinates of the input point (y1, z1) and root_flag = True (if point is found).
# R = 1.
def func(y2, z2, m, init_approx):
    y1 = 0
    z1 = 0
    count = 0
    root_flag = False

    sin = init_approx
    cos = (1 - sin ** 2) ** 0.5
    tmp1 = (m ** 2 - sin ** 2) ** 0.5
    g = sin * (tmp1 - cos) / (sin ** 2 + cos * tmp1)

    while count < max_iterations:
        # According to the estimate of the right function (g) find a new sine (sin_new):
        tmp2 = y2 + z2 * g
        tmp3 = 1 + g ** 2
        sin_new = (tmp2 + g * (tmp3 - tmp2 ** 2) ** 0.5) / tmp3
        # Optimization for the case when the algorithm diverges:
        if sin_new > 1 or sin_new < 0:
            break

        # I find what for a given sine (sin_new) the function on the right is equal:
        cos_new = (1 - sin_new ** 2) ** 0.5
        tmp1_new = (m ** 2 - sin_new ** 2) ** 0.5
        g_new = sin_new * (tmp1_new - cos_new) / (sin_new ** 2 + cos_new * tmp1_new)

        # I find what for a given sine (sin_new) the function on the left is equal:
        fnew = (sin_new - y2) / (z2 + cos_new)
        d = abs(fnew - g_new)
        if d < tolerance:
            y1 = sin_new
            root_flag = True
            break
        count += 1
        g = g_new

    if root_flag and y2 <= y1 <= 1:
        z1 = - (1 - y1 ** 2) ** 0.5
        if z2 < z1 or 0 < z1:
            root_flag = False
    else:
        root_flag = False

    return root_flag, y1, z1


# The function solves the equation and finds the right root.
# R = 1.
def func2(y2, z2, m, init_approx):
    y1 = 0
    z1 = 0
    count = 0
    root_flag = False

    sin = init_approx
    cos = (1 - sin ** 2) ** 0.5
    f = (sin - y2) / (z2 + cos)

    while count < max_iterations:
        # According to the assessment of the left function (f) find a new sine (sin_new):
        tmp1 = f ** 2 + 1
        tmp2 = tmp1 ** 0.5
        sin_new = f * m / ((m ** 2 + 1) * tmp1 - 2 * m * tmp2) ** 0.5
        # Optimization for the case when the algorithm diverges:
        if sin_new > 1 or sin_new < 0:
            break

        # I find what for a given sine (sin_new) the function on the left is equal:
        cos_new = (1 - sin_new ** 2) ** 0.5
        f_new = (sin_new - y2) / (z2 + cos_new)

        # I find what for a given sine (sin_new) the function on the right is equal:
        tmp3 = (m ** 2 - sin_new ** 2) ** 0.5
        f2 = sin_new * (tmp3 - cos_new) / (sin_new ** 2 + cos_new * tmp3)
        d = abs(f2 - f_new)
        if d < tolerance:
            y1 = sin_new
            root_flag = True
            break
        count += 1
        f = f_new

    if root_flag and y2 <= y1 <= 1:
        z1 = - (1 - y1 ** 2) ** 0.5
        if z2 < z1 or 0 < z1:
            root_flag = False
    else:
        root_flag = False

    return root_flag, y1, z1


def iterative_method(y2, z2, m, cur_region):
    info = "not_success"
    y1 = 0
    z1 = 0
    y1_2 = 0
    z1_2 = 0
    if cur_region == "one_root":
        root_flag, y1, z1 = func(y2, z2, m, y2)
        if root_flag:
            if 0 <= y1 <= 1 and -1 <= z1 <= 0:
                info = "success"
            else:
                print("Error in iterative_method()! 0 <= y1 <= 1 and -1 <= z1 <= 1 not performed.")
                info = "not_success"
    elif cur_region == "two_roots":
        root_flag, y1, z1 = func(y2, z2, m, y2)
        root_flag2, y1_2, z1_2 = func2(y2, z2, m, 1)
        if root_flag and root_flag2:
            if 0 <= y1 <= y1_2 <= 1 and -1 <= z1 <= z1_2 <= 0:
                info = "success"
            else:
                print(
                    "Error in iterative_method()! 0 <= y1 <= y1_2 <= 1 and -1 <= z1 <= z1_2 <= 0 not performed.")
                info = "not_success"
    elif cur_region == "no_root":
        info = "not_success"
    else:
        print("Error iterative_method() function!")
    return info, y1, z1, y1_2, z1_2