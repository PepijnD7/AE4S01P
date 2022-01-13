import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize
from Read_Provided_data import plot_given_data

# TODO: Implement variation of a based on temperature variations


# set given to True for validation
given = False

# Propellant data
n = 0.222
a = 0.005132 * (10 ** (-6)) ** n
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

const = [d_port, d_out, l_p, alphast, eps, a, n, P_a, m_p, rho_p]
conII = [d_port, d_out, l_p, alphaII, eps, a, n, P_a, m_p, rho_p]
const_given = [d_port_given, d_out_given, l_p_given, alpha_given, eps, a, n, P_a, m_p_given, rho_p_given]


# Regression Rate
def regrate(P, a, n):
    r = a * P**n
    return r


def DivLoss(alpha):
    return (1 + np.cos(alpha)) / 2


def Kerckhove(Gamma):
    return np.sqrt(Gamma) * (2 / (1 + Gamma))**(0.5 * (Gamma + 1) / (Gamma - 1))

def linearizeAccumulate(param_str, p):
    data_RPA = {'Pressure': [0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1],
                'Gas Constant' : [217.5 , 214.9 , 213.5 , 212.5 , 211.7 , 211.2 , 210.7 , 210.3 , 209.9 , 209.6],
                'Temperature': [1451.92 , 1491.25 , 1512.33 , 1528.3707 , 1539.5906 , 1548.4226 , 1555.6289 , 1561.6622 , 1566.8127 , 1571.2774],
                'Gamma': [1.1767 , 1.1658 , 1.1598 , 1.1558 , 1.1529 , 1.1506 , 1.1488 , 1.1472 , 1.1460 , 1.1449],
                'Characteristic velocity_opt': [892.79 , 898.35 , 901.27 , 903.17 , 904.53 , 905.56 , 906.39 , 907.06 , 907.62 , 908.1],
                'Thrust coefficient_opt':[1.4237 , 1.4241 , 1.4244 , 1.4246 , 1.4246 , 1.4247 , 1.4248 , 1.4248 , 1.4248 , 1.4248]
                }

    p /= 1000000
    p_lst = data_RPA['Pressure']

    i = 0
    while p_lst[i] < p:
        i += 1

    dif = p - p_lst[i-1]
    step = p_lst[i] - p_lst[i-1]

    smaller = data_RPA[param_str][i-1]
    bigger  = data_RPA[param_str][i]

    return smaller + (bigger - smaller) * (dif / step)


def FindPratio(Gamma, eps):
    Guess = 0.1
    for _ in range(10):
        Guess = (Kerckhove(Gamma)**2 / (eps**2 * (2 * Gamma / (Gamma - 1)) * (1 - Guess**((Gamma - 1) / Gamma))))**(Gamma / 2)
    return Guess

dt = 0.005
# Simulation
def Simulation(con):
    g0 = 9.81
    d_port, d_out, l_p, alpha, eps, a, n, P_a, m_p, rho_p = con
    a = (1.55 * 10 ** (-5) * (T_a - 273.15) + 4.47 * 10 ** (-3))* (10 ** (-6)) ** n
    # INITIAL CONDITIONS and CHAMBER FILL:

    # Parameters that do not change:
    A_t = (d_t ** 2 * np.pi / 4)
    V_p = l_p * ((d_out/2)**2 - (d_port/2)**2) * np.pi
    rho_p = m_p/V_p
    dt = 0.005

    # Parameters that will change:
    V_c = l_p * (d_port/2)**2 * np.pi
    P_c = 1*10**5
    S = np.pi * l_p * d_port
    dpdt = 1
    m_list = []
    Isp_list = []
    p_list = []
    T_list = []
    r_list = []
    I_list = [0]
    pepa_list = []

    m_out = 0
    m_in = regrate(P_c, a, n) * rho_p * S
    accumdiff = 1

    # while np.abs(accumdiff)>0.01:
    # while np.abs(accumdiff)>0.01:
    #     print("in sim pc1: ",P_c)
    #     S = np.pi * l_p * d_port
    #     Gamma = linearizeAccumulate('Gamma', P_c)
    #     vdk = Kerckhove(Gamma)
    #     c_star = linearizeAccumulate('Characteristic velocity_opt', P_c)
    #     C_f0 = linearizeAccumulate('Thrust coefficient_opt', P_c)
    #     T_c = linearizeAccumulate("Temperature", P_c)
    #     R = linearizeAccumulate('Gas Constant', P_c)
    #     print("in sim pc2: ",P_c)
    #
    #     rho_c = P_c/R/T_c
    #     r = regrate(P_c, a, n)
    #     m_in = r * S
    #     m_out = A_t * P_c / c_star
    #     accumdiff = m_out - m_in*(rho_p - rho_c)
    #     print("ACUUUUUUM",np.abs(accumdiff), m_in*(rho_p - rho_c), m_out)
    #     dpdt = vdk**2/V_c * (c_star**2 * (rho_p - rho_c) * m_in - m_out)  # Add gas density
    #     C_f = C_f0 * DivLoss(alpha) + (FindPratio(Gamma, eps) - P_a / P_c) * eps
    #     T = C_f * P_c * A_t
    #     Isp = T / g0 / m_out
    #     m_list.append(m_out)
    #     Isp_list.append(Isp)
    #     p_list.append(P_c)
    #     pepa_list.append(P_c * FindPratio(Gamma, eps) / P_a)
    #     T_list.append(T)
    #     r_list.append(r)
    #     prev_I = I_list[-1]
    #     I_list.append(prev_I + (T * dt))
    #
    #     d_port+= r*dt
    #     P_c += dpdt*dt
    #     print("DPDT",dpdt)

    for _ in range(10):
        c_star = linearize('Characteristic velocity_opt', P_c)
        P_c = (c_star * rho_p * a * S / A_t)**(1 / (1-n))
    r = regrate(P_c, a, n)
    m_dot = r * S * rho_p

    c_star = linearize('Characteristic velocity_opt', P_c)
    Gamma = linearize('Gamma', P_c)
    C_f0 = linearize('Thrust coefficient_opt', P_c)
    T_c = linearize("Temperature", P_c)
    R = linearize('Gas constant', P_c) * 1000
    rho_c = P_c / R / T_c

    C_f = C_f0 * DivLoss(alpha) + (FindPratio(Gamma, eps) + P_a / P_c) * eps
    o = 0
    # BURN LOOP

    while d_port <= d_out:
        o += 1
        #print("Running", o, P_c)
        S_burn = d_port * np.pi * l_p
        P_c = (c_star * (rho_p-rho_c) * a * S_burn / A_t) ** (1 / (1 - n))
        T_c = linearize("Temperature", P_c)
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

    return np.array(t_list), np.array(p_list), np.array(I_list), np.array(m_list), np.array(T_list), \
           np.array(r_list), np.array(Isp_list), np.array(pepa_list)


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

    # plt.plot(t_II, pepa_II, label='Config. II')
    # plt.xlabel("Time [s]", fontsize=13)
    # plt.ylabel(r"$p_e/p_a$", fontsize=13)
    # plt.grid()
    # plt.show()


if given:
    plt.plot(t_given, p_given * 10 ** -6, label='Chamber pressure')
    plt.xlabel("Time [s]", fontsize=13)
    plt.ylabel("Chamber pressure [MPa]", fontsize=13)
    plt.title("Chamber pressure", fontsize=13)
    plt.grid()
    plt.show()

    plot_given_data([t_given, p_given])
