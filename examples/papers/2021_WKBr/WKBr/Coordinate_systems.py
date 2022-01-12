from math import atan2, sin, cos, asin


# Function for finding the coordinates of a point in the meridian plane.
# M.p. - this is the plane formed by the entry point and the z-axis.
def coordinates_in_meridional_plane(x_lab, y_lab, z_lab):
    return (x_lab ** 2 + y_lab ** 2) ** 0.5, z_lab


# Function to find the angle of rotation between the meridian plane and the laboratory coordinate system.
# The angle is counted from the OX axis counterclockwise.
# Angle range [-pi, pi].
def find_rotation_angle_in_mp(x_lab, y_lab):
    return atan2(y_lab, x_lab)


# Function for finding point coordinates in lab.c.s.
def coordinates_in_lab_plane(y_mp, z_mp, rotation_angle):
    x_lab = y_mp * cos(rotation_angle)
    y_lab = y_mp * sin(rotation_angle)
    z_lab = z_mp
    return x_lab, y_lab, z_lab


# Function for finding the angle i-t in lab.c.s. in section xz:
def delta_angle(x1, m):
    i = asin(x1)
    t = asin(x1 / m)
    return i - t


# Function for finding the cosine of the angle of passage inside the sphere in the meridian plane:
def cos_refracted_angle_mp(y1, m):
    sin_ = y1 / m
    cos_ = (1 - sin_ ** 2) ** 0.5
    return cos_


# Function for rotating the electric field vector due to the curvature of the particle surface.
def apply_rotation_electric_field_vector(exr, exi, da):
    exr_new = exr * cos(da)
    exi_new = exi * cos(da)
    ezr_new = exr * sin(da)
    ezi_new = exi * sin(da)
    return exr_new, exi_new, 0, 0, ezr_new, ezi_new


# The sum of electric fields in the region of two solutions.
def sum_ef(exr1, exr2, exi1, exi2,
           eyr1, eyr2, eyi1, eyi2,
           ezr1, ezr2, ezi1, ezi2):
    exr_new = exr1 + exr2  # the sum of the real X-component of the electric field
    exi_new = exi1 + exi2  # the sum of the imaginary X-component of the electric field
    eyr_new = eyr1 + eyr2
    eyi_new = eyi1 + eyi2
    ezr_new = ezr1 + ezr2
    ezi_new = ezi1 + ezi2
    return exr_new, exi_new, eyr_new, eyi_new, ezr_new, ezi_new


# Function for finding the electric field vector Epar, Eper.
def electric_field_vectors_in_mp(exr, exi, rotation_angle):
    er_par = exr * cos(rotation_angle)
    ei_par = exi * cos(rotation_angle)
    er_per = exr * sin(rotation_angle)
    ei_per = exi * sin(rotation_angle)
    return er_par, ei_par, er_per, ei_per


# Function for finding the electric field vector Exr, Exi, Eyr, Eyi.
def electric_field_vectors_in_lab(er_par, ei_par, er_per, ei_per, rotation_angle):
    exr = er_par * cos(rotation_angle) + er_per * sin(rotation_angle)
    exi = ei_par * cos(rotation_angle) + ei_per * sin(rotation_angle)
    eyr = 0  # - er_par * sin(rotation_angle) + er_per * cos(rotation_angle)
    eyi = 0  # - ei_par * sin(rotation_angle) + ei_per * cos(rotation_angle)
    return exr, exi, eyr, eyi


# Function for applying the transmission coefficient.
def apply_transmission_coefficient(exr, exi, t_per, t_par, rotation_angle):
    er_par, ei_par, er_per, ei_per = electric_field_vectors_in_mp(exr, exi, rotation_angle)
    er_parc = complex(er_par, ei_par)
    er_perc = complex(er_per, ei_per)
    er_parc *= t_par
    er_perc *= t_per
    exr_new, exi_new, eyr_new, eyi_new = electric_field_vectors_in_lab(er_parc.real, er_parc.imag,
                                                                       er_perc.real, er_perc.imag,
                                                                       rotation_angle)
    return exr_new, exi_new, eyr_new, eyi_new
