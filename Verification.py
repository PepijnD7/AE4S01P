import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize

# TODO: Implement variation of a based on temperature variations

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

const = [d_port, d_out, l_p, alphast, eps, a, n, P_a, Pref, m_p, rho_p]
conII = [d_port, d_out, l_p, alphaII, eps, a, n, P_a, Pref, m_p, rho_p]


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

# plt.plot(P_list, r_list*10**3)
# plt.xlabel("Chamber Pressure (Pa)")
# plt.ylabel("Regression rate (mm/s)")
# plt.grid()
# plt.show()

r_c = r_list[P_list == 100*10**6]
t_list = np.arange(0, 10, 0.1)
d_list, S_list, V_list, m_list = burnsurf(r_c, l_p, d_port, d_out, t_list)
d_reg = d_list-d_port

# fig, ax1 = plt.subplots()
#
# ax1.set_xlabel("Distance regressed (mm)")
# ax1.set_ylabel("Propellant mass (kg)")
# ax1.grid()
# ax1.plot(d_reg*10**3, m_list, color='r', label='Propellant mass')
#
# ax2 = ax1.twinx()
#
# ax2.set_xlabel("Distance regressed (mm)")
# ax2.set_ylabel("Burn surface (m^2)")
# ax2.plot(d_reg*10**3, S_list, color='b', label='Burn Surface')
# fig.tight_layout()
# fig.legend()
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

    a = 0.005132 * np.exp(olde * (T_a - 273.15))
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
    first = True
    while d_port <= d_out:
        S_burn = d_port * np.pi * l_p

        a *= (10 ** (-6)) ** n
        P_c = (c_star * rho_p * a * S_burn / A_t) ** (1 / (1 - n))
        T_c = linearize("Temperature", P_c)

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
        if first:
            print("INITIAL CONDITIONS LOOP 2")
            print("c_star=", c_star)
            print("d_p=", d_port)
            print("rho_c=", c_star)
            print("I=", I_list)
            first = False

    t_list = np.arange(0, dt*len(p_list), dt)
    I_list = I_list[1:]

    return np.array(t_list), np.array(p_list), np.array(I_list), np.array(m_list), np.array(T_list), np.array(r_list), np.array(Isp_list), np.array(pepa_list)


t_II, p_II, I_II, m_II, T_II, r_II, Isp_II, pepa_II = Simulation(conII)

p_hand, m_hand, r_hand, T_hand, Isp_hand, I_hand = 1367496.858, 0.082722, 0.00550162, 99.91351429, 123.1214536, 19.18135
# Printing output


print("CHAMBER PRESSURE:")
print("Code:", p_II[1], "\n")
print("Difference:", (p_hand-p_II[1])/p_II[1]*100, "\n")


print("MASS FLOW:")
print("Code:", m_II[1], "\n")
print("Difference:", (m_hand-m_II[1])/m_II[1]*100, "\n")

print("REGRESSION RATE:")
print("Code:", r_II[1], "\n")
print("Difference:", (r_hand-r_II[1])/r_II[1]*100, "\n")

print("THRUST:")
print("Code:", T_II[1], "\n")
print("Difference:", (T_hand-T_II[1])/T_II[1]*100, "\n")

print("SPECIFIC IMPULSE:")
print("Code:", Isp_II[1], "\n")
print("Difference:", (Isp_hand-Isp_II[1])/Isp_II[1]*100, "\n")

print("TOTAL IMPULSE:")
print("Code:", I_II[1], "\n")
print("Difference:", (I_hand-I_II[1])/I_II[1]*100, "\n")
