from math import pi, e, sqrt, tan, log

print("All values and formulas are in SI units except for spectral line wavelengths, which are in angstroms.")

# UNIT CONVERSIONS
ly_m = 9.461 * 1e15
pc_m = 3.086 * 1e16
au_m = 1.496 * 1e11
m_ly = 1 / ly_m
m_pc = 1 / pc_m
m_au = 1 / au_m
pc_ly = pc_m / ly_m
pc_au = pc_m / au_m
ly_pc = 1 / pc_ly
au_pc = 1 / pc_au
ly_au = ly_m / au_m
au_ly = 1 / ly_au

J_eV = 6.242 * 1e18
ev_J = 1 / J_eV

# IMPORTANT NUMBERS AND CONSTANTS (SI UNITS)

G = 6.67 * 1e-11

u_0 = (4 * pi) * 1e-7
e_0 = 8.85418 * 1e-12
c = 1 / sqrt(u_0 * e_0)

hbar = 1.05457*1e-34
h = hbar * 2 * pi

a_0 = 5.292 * 1e-11
alpha = 1 / 137

u = 1.66 * 1e-27
q_e = 1.602 * 1e-19
m_e = 9.109 * 1e-31
m_p = 1.673 * 1e-27
m_n = 1.675 * 1e-27

N_A = 6.022 * 1e23
k_B = 1.381 * 1e-23
R = N_A * k_B

sigma = 5.67037 * 1e-8
b = 2.8977 * 1e-3

H_0 = 2.25 * 1e-18

R_H = 1.097 * 1e7

### DATA ###
M_sun = 1.989 * 1e30
R_sun = 6.96 * 1e8
L_sun = 3.9 * 1e26
T_sun = 5777


emission_lines = {"O VI": 1033.82, "N V": 1240.81, "O I": [1305.53, 6302.05, 6365.54, ], "C II": [1335.31, 2326.0], "Si IV": 1397.61, "Si IV + O IV": 1399.8, "C IV": 1549.48, "He II": 1640.4, "O III": [
    1665.85, 4364.44, 4932.6, 4960.3, 5008.24], "Al III": 1857.4, "C III": 1908.73, "Ne IV": 2439.5, "Mg II": 2799.12, "Ne V": 3346.79, "O II": [3727.09, 3729.86], "He I": 3889.0, "S II": [4072.3, 6718.29, 6732.67], "N I": 6529.03, "N II": 6585.27}

absorption_lines = {"K": 3934.78, "H": 3969.59, "G": 4305.61,
                    "Mg": 5176.7, "Na": 5895.6, "Ca II": [8500.36, 8544.44, 8664.52]}


############# FUNCTIONS ##############
def quad(a, b, c):
    d = b**2 - 4 * a * c
    return (-b + sqrt(d)) / 2 / a, (-b - sqrt(d)) / 2 / a


def orbital_speed(M, r, a):
    return sqrt(G * M * (2.0 / r - 1.0 / a))


def escape_speed(M, r):
    return sqrt(2 * G * M / r)


def swarzchild_radius(M):
    return 2 * G * M / c ** 2


def kepler3law(M, P=None, a=None):
    if P:
        return (P ** 2 * G * M / 4 / pi ** 2) ** (1/3)
    else:
        return (4 * pi ** 2 * a ** 3 / G / M) ** 0.5


def harm_mean(a, b):
    return a * b / (a + b)


def redshift(original, final):
    return (final - original)/original


def doppler_shift(beta):
    return sqrt((1 + beta)/(1 - beta))


def redshift_vel(original=None, final=None, z=None):
    if not z:
        z = redshift(original, final)
    return c * ((z + 1) ** 2 - 1)/((z + 1) ** 2 + 1)


def hubbles_law(original=None, final=None, z=None, v=None):
    if original and final:
        z = redshift(original, final)
        v = redshift_vel(z=z)
    elif z:
        v = redshift_vel(z=z)

    return v / H_0


def parallax(p, a=None, SI=False):
    if SI:
        return a / tan(p)
    else:
        print("Answer in PARSECS (PC), input should have been in arcseconds ('')")
        return 1 / p


def planck_func(T, lambda_=None, nu=None):
    if lambda_:
        return 2 * h * c ** 2 / lambda_ ** 5 / (e ** (h * c / lambda_ / k_B / T) - 1)
    else:
        return 2 * h * nu ** 3 / c ** 2 / (e ** (h * nu / k_B / T) - 1)


def hydrogen_lines(m, n):
    return 1.0 / R_H / (1.0 / m ** 2.0 - 1.0 / n ** 2.0)


def distance(distance_mod):
    print("Output is in PARSECS")
    return 10.0 * 10.0 ** ((distance_mod) / 5.0)


def saha_eq(g_i, g_i1, T, E_i, n_e=None, P_e=None):
    if n_e:
        return 2 * g_i1/g_i/n_e * (2 * pi * m_e * k_B * T / h**2)**1.5 * e**(-E_i / k_B / T)
    else:
        return 2 * k_B * T * g_i1/g_i/P_e * (2 * pi * m_e * k_B * T / h**2)**1.5 * e**(-E_i / k_B / T)


def cepheid(T):
    return -2.76 * (log(T)/log(10) - 1) - 4.16

# TODO: Jean's Mass, Jeans Radiuss
