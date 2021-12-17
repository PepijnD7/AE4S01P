import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize
from Read_Provided_data import plot_given_data

# TODO: Implement variation of a based on temperature variations


# set given to True for validation
given = False

# Propellant data
a = 0.005132
n = 0.222
Pref = 1*10**6
P_a = 101325
rho_p = 1720.49
m_p = 0.758
l_p = 0.107
d_out = 0.0766
d_port = 0.025
d_t = 8.37 / 1000
alphast = 12 * np.pi / 180          # Standard setup
alphaII = 15 * np.pi / 180          # Configuration II
eps = 4
T_a = 273.15 + 15

# Given data
l_p_given = 99.60 / 1000
d_out_given = 76.6 / 1000
d_port_given = 24.90 / 1000
m_p_given = 741.8 / 1000
alpha_given = 12 * np.pi / 180
rho_p_given = 1786.5

const = [d_port, d_out, l_p, alphast, eps, a, n, P_a, Pref, m_p, rho_p]
conII = [d_port, d_out, l_p, alphaII, eps, a, n, P_a, Pref, m_p, rho_p]
const_given = [d_port_given, d_out_given, l_p_given, alpha_given, eps, a, n, P_a, Pref, m_p_given, rho_p_given]


# Regression Rate
def regrate(P, a, n):
    r = a * (P/Pref)**n
    return r


def DivLoss(alpha):
    return (1 + np.cos(alpha)) / 2


def Kerckhove(Gamma):
    return np.sqrt(Gamma) * (2 / (1 + Gamma))**(0.5 * (Gamma + 1) / (Gamma - 1))


def FindPratio(Gamma, eps):
    Guess = 0.1
    for _ in range(10):
        Guess = (Kerckhove(Gamma)**2 / (eps**2 * (2 * Gamma / (Gamma - 1)) * (1 - Guess**((Gamma - 1) / Gamma))))**(Gamma / 2)
    return Guess


P_list = np.arange(1, 500)*10**6
r_list = regrate(P_list, a, n)


# Propellant mass variation
V_p = l_p * ((d_out/2)**2 - (d_port/2)**2) * np.pi
rho_p = m_p/V_p


def burnsurf(r, l, d_i, d_o, t):
    d_ilist = np.ones((len(t)))*d_i
    d_list = d_ilist + 2*r*t
    d_list[d_list>=d_o] = d_o
    S_list = np.pi * l * d_list
    V_list = l * ((d_o/2)**2 - (d_list/2)**2) * np.pi
    m_list = V_list * rho_p
    return d_list, S_list, V_list, m_list


# Simulation

# plt.subplot(1, 3, 1)
# plt.plot(P_list, r_list*10**3, label='Regression rate')
# plt.xlabel("Chamber Pressure (Pa)")
# plt.ylabel("Regression rate (mm/s)")
# plt.grid()
# plt.legend()
# plt.show()

r_c = r_list[P_list == 100*10**6]
t_list = np.arange(0, 10, 0.1)
d_list, S_list, V_list, m_list = burnsurf(r_c, l_p, d_port, d_out, t_list)
d_reg = d_list-d_port

# fig, ax1 = plt.subplots()
# plt.subplot(1, 3, 2)
# plt.xlabel("Distance regressed (mm)")
# plt.ylabel("Propellant mass (kg)")
# plt.grid()
# plt.plot(d_reg*10**3, m_list, label='Propellant mass')
# plt.legend()

# ax2 = ax1.twinx()

# plt.subplot(1, 3, 3)
# plt.xlabel("Distance regressed (mm)")
# plt.ylabel("Burn surface (m^2)")
# plt.plot(d_reg*10**3, S_list, label='Burn surface')
# # fig.tight_layout()
# plt.legend()
# plt.grid()
# plt.show()


# Initial conditions
# P_c = 5 * Pref
# S = np.pi * l_p * d_port
# A_t = (d_t ** 2 * np.pi / 4)
# a *= (10 ** (-6)) ** n
#
# for _ in range(10):
#     c_star = linearize('Characteristic velocity_opt', P_c)
#     # print('c_star: ', c_star)
#
#     P_c = (c_star * rho_P * a * S / A_t)**(1 / (1-n))
#     # print('P_c: ', P_c)
#
# a = 0.005132
# r = regrate(P_c, a, n)
# m_dot = r * S * rho_P
#
# c_star = linearize('Characteristic velocity_opt', P_c)
# Gamma = linearize('Gamma',P_c)
# C_f0 = linearize('Thrust coefficient_opt', P_c)
# C_f = C_f0 * DivLoss(alpha) + (FindPratio(Gamma,eps) + P_a / P_c) * eps
# T = C_f * P_c * A_t

# print(P_c, r, m_dot, T, Gamma, C_f0, C_f)


# Simulation
def Simulation(con):
    g0 = 9.81
    olde = 0
    d_port, d_out, l_p, alpha, eps, a, n, P_a, Pref, m_p, T_a = con

    V_p = l_p * ((d_out/2)**2 - (d_port/2)**2) * np.pi
    rho_p = m_p/V_p

    # Initial conditions
    P_c = 5 * Pref
    S = np.pi * l_p * d_port
    A_t = (d_t ** 2 * np.pi / 4)
    a *= (10 ** (-6)) ** n

    for _ in range(10):
        c_star = linearize('Characteristic velocity_opt', P_c)
        P_c = (c_star * rho_p * a * S / A_t)**(1 / (1-n))

    a = 0.005132
    r = regrate(P_c, a, n)
    m_dot = r * S * rho_p

    c_star = linearize('Characteristic velocity_opt', P_c)
    Gamma = linearize('Gamma', P_c)
    C_f0 = linearize('Thrust coefficient_opt', P_c)
    C_f = C_f0 * DivLoss(alpha) + (FindPratio(Gamma, eps) + P_a / P_c) * eps
    T = C_f * P_c * A_t

    # BURN LOOP
    dt = 0.1

    m_list = []
    Isp_list = []
    p_list = []
    T_list = []
    r_list = []
    I_list = [0]
    pepa_list = []

    while d_port <= d_out:
        S_burn = d_port * np.pi * l_p

        a *= (10 ** (-6)) ** n
        P_c = (c_star * rho_p * a * S_burn / A_t) ** (1 / (1 - n))
        T_c = linearize("Temperature",P_c)

        a = 0.005132 * np.exp(olde * (T_a - 273.15))
        # a = (1.55*10**(-5) * (T_a - 273.15) + 4.47*10**(-3))
        # print(T_c,a)
        r_rate = regrate(P_c, a, n)

        m_dot = rho_p * r_rate * S_burn
        c_star = linearize('Characteristic velocity_opt', P_c)
        Gamma = linearize('Gamma', P_c)
        C_f0 = linearize('Thrust coefficient_opt', P_c)
        C_f = C_f0 * DivLoss(alpha) + (FindPratio(Gamma, eps) - P_a / P_c) * eps
        T = C_f * P_c * A_t
        Isp = T / g0 / m_dot
        # print(Isp - linearize("Specific impulse (by mass)", P_c))
        m_list.append(m_dot)
        Isp_list.append(Isp)
        p_list.append(P_c)
        pepa_list.append(P_c*FindPratio(Gamma, eps)/P_a)
        T_list.append(T)
        r_list.append(r_rate)
        prev_I = I_list[-1]
        I_list.append(prev_I + (T * dt))
        d_port += 2 * r_rate * dt

    t_list = np.arange(0, dt*len(p_list), dt)
    I_list = I_list[1:]

    return np.array(t_list), np.array(p_list), np.array(I_list), np.array(m_list), np.array(T_list), np.array(r_list), np.array(Isp_list), np.array(pepa_list)


t_st, p_st, I_st, m_st, T_st, r_st, Isp_st, pepa_st = Simulation(const)
t_II, p_II, I_II, m_II, T_II, r_II, Isp_II, pepa_II = Simulation(conII)
t_given, p_given, _, _, _, _, _, _ = Simulation(const_given)


if not given:
    # Printing output
    print("CHAMBER PRESSURE:\n")
    print("Max [MPa]:", np.max(p_II)*10**-6)
    print("Average [MPa]:", np.average(p_II)*10**-6, "\n")

    print("MASS FLOW:\n")
    print("Max:", np.max(m_II))
    print("Ave:", np.average(m_II), "\n")

    print("REGRESSION RATE:\n")
    print("Max [mm/s]:", np.max(r_II)*10**3)
    print("Ave [mm/s]:", np.average(r_II)*10**3, "\n")

    print("THRUST:\n")
    print("Max:\n 12 deg:", np.max(T_st), "\n 15 deg:", np.max(T_II), "\n Compare [%]:",
          (np.max(T_II)-np.max(T_st))/np.max(T_st)*100, "\n")

    print("Ave:\n 12 deg:", np.average(T_st), "\n 15 deg:", np.average(T_II), "\n Compare [%]:",
          (np.average(T_II)-np.average(T_st))/np.average(T_st)*100, "\n")

    print("SPECIFIC IMPULSE:\n")
    print("Max:\n 12 deg:", np.max(Isp_st), "\n 15 deg:", np.max(Isp_II), "\n Compare [%]:",
          (np.max(Isp_II)-np.max(Isp_st))/np.max(Isp_st)*100, "\n")

    print("Ave:\n 12 deg:", np.average(Isp_st), "\n 15 deg:", np.average(Isp_II), "\n Compare [%]:",
          (np.average(Isp_II)-np.average(Isp_st))/np.average(Isp_st)*100, "\n")

    print("TOTAL IMPULSE:\n")
    print("12 deg:", np.max(I_st), "\n15 deg:", np.max(I_II), "\n Compare [%]:",
          (np.max(I_II)-np.max(I_st))/np.max(I_st)*100, "\n")


    # Plotting outputs
    plt.subplot(1, 3, 1)
    plt.plot(t_st, p_II*10**-6, label='Chamber pressure')
    plt.xlabel("Time [s]", fontsize=13)
    plt.ylabel("Chamber pressure [MPa]", fontsize=13)
    plt.title("Chamber pressure", fontsize=13)
    plt.grid()

    plt.subplot(1, 3, 2)
    plt.plot(t_II, m_II, label='Mass flow')
    plt.xlabel("Time [s]", fontsize=13)
    plt.ylabel("Mass flow [kg/s]", fontsize=13)
    plt.title("Mass flow", fontsize=13)
    plt.grid()

    plt.subplot(1, 3, 3)
    plt.plot(t_II, r_II*10**3, label='Regression rate')
    plt.xlabel("Time [s]", fontsize=13)
    plt.ylabel("Regression rate [mm/s]", fontsize=13)
    plt.title("Regression rate", fontsize=13)
    plt.grid()
    plt.show()


    plt.subplot(1, 3, 1)
    plt.plot(t_II, T_II, label='Config. II')
    plt.plot(t_st, T_st, label='Standard')
    plt.xlabel("Time [s]", fontsize=13)
    plt.ylabel("Thrust [N]", fontsize=13)
    plt.title("Thrust", fontsize=13)
    plt.grid()
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(t_II, Isp_II, label='Config. II')
    plt.plot(t_st, Isp_st, label='Standard')
    plt.xlabel("Time [s]", fontsize=13)
    plt.ylabel("Specific impulse [s]", fontsize=13)
    plt.title("Specific impulse", fontsize=13)
    plt.grid()
    plt.legend()

    plt.subplot(1, 3, 3)
    plt.plot(t_II, I_II, label='Config. II')
    plt.plot(t_st, I_st, label='Standard')
    plt.xlabel("Time [s]", fontsize=13)
    plt.ylabel("Total impulse [Ns]", fontsize=13)
    plt.title("Total impulse", fontsize=13)
    plt.grid()
    plt.legend()
    plt.show()

    plt.plot(t_II, pepa_II, label='Config. II')
    plt.xlabel("Time [s]", fontsize=13)
    plt.ylabel(r"$p_e/p_a$", fontsize=13)
    plt.grid()
    plt.show()


if given:
    plt.plot(t_given, p_given * 10 ** -6, label='Chamber pressure')
    plt.xlabel("Time [s]", fontsize=13)
    plt.ylabel("Chamber pressure [MPa]", fontsize=13)
    plt.title("Chamber pressure", fontsize=13)
    plt.grid()
    plt.show()

    plot_given_data([t_given, p_given])
