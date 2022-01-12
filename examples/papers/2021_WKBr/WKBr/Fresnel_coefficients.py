

# The function finds the transmission coefficient for the case of an X-polarized wave,
# calculations in the meridial plane (y, z).
# I use the function when there is an imaginary refractive index (attenuation in the medium).
# Formulas from Peter C.Y. Chang et al "Ray tracing in absorbing media" (2005).
def transmission_coefficient(y1, mr, mi, k0):
    sini = y1
    sint = y1 / mr
    cosi = (1 - sini ** 2) ** 0.5
    cost = (1 - sint ** 2) ** 0.5
    if mi == 0:
        tper = 2 * cosi / (cosi + mr * cost)
        tpar = 2 * cosi / (cost + mr * cosi)
        if tper < 0 or tper > 1:
            print("Error in transmission_coefficient() function! t_per is invalid.")
        if tpar < 0 or tpar > 1:
            print("Error in transmission_coefficient() function! t_par is invalid.")
    elif mi > 0:
        ki = k0 * cosi
        sinti = 0 # sini / mi
        costi = 1 # (1 - sinti ** 2) ** 0.5
        kt = k0 * complex(mr * cost, mi * costi)
        tper = 2 * ki / (ki + kt)
        n1 = 1
        n2 = complex(mr, mi)
        nratio = n2 / n1
        nratio2 = nratio ** 2
        tpar = 2 * ki * nratio / (nratio2 * ki + kt)
    else:
        print("Error in optical_len() function! mi is not valid.")
    return tper, tpar


# Function for finding effective real (N) and imaginary (K) refractive index:
# mr - real,
# mi - imaginary.
# Formulas from Peter C.Y. Chang et al. "Ray tracing in absorbing media" (2005).
def effective_refractive_indices(mr, mi, cos_t):
    mr2 = mr ** 2
    mi2 = mi ** 2
    N2 = (mr2 - mi2 + ((mr2 - mi2) ** 2 + 4 * mr2 * mi2 / (cos_t ** 2)) ** 0.5) / 2
    N = N2 ** 0.5
    K = (N2 - mr2 + mi2) ** 0.5
    return N, K


