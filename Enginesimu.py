import numpy as np
from matplotlib import pyplot as plt

# Propellant data
a = 0.005132
n = 0.222
Pref = 1*10**6
m_p = 0.758
l_p = 0.107
d_out = 0.0766
d_port = 0.025

# Regression Rate


def regrate(P, a, n):
    r = a * (P/Pref)**n
    return r


P_list = np.arange(1, 500)*10**6
r_list = regrate(P_list, a, n)


# Propellant mass variation
V_p = l_p * ((d_out/2)**2 - (d_port/2)**2) * np.pi
rho_P = m_p/V_p


def burnsurf(r, l, d_i, d_o, t):
    d_ilist = np.ones((len(t)))*d_i
    d_list = d_ilist + r*t
    d_list[d_list>=d_o] = d_o
    S_list = np.pi * l * d_list
    V_list = l * ((d_o/2)**2 - (d_list/2)**2) * np.pi
    m_list = V_list * rho_P
    return d_list, S_list, V_list, m_list


# Simulation

plt.plot(P_list, r_list*10**3)
plt.xlabel("Chamber Pressure (Pa)")
plt.ylabel("Regression rate (mm/s)")
plt.grid()
plt.show()

r_c = r_list[P_list == 100*10**6]
t_list = np.arange(0, 10, 0.1)
d_list, S_list, V_list, m_list = burnsurf(r_c, l_p, d_port, d_out, t_list)
d_reg = d_list-d_port

fig, ax1 = plt.subplots()

ax1.set_xlabel("Distance regressed (mm)")
ax1.set_ylabel("Propellant mass (kg)")
ax1.grid()
ax1.plot(d_reg*10**3, m_list, color='r', label='Propellant mass')

ax2 = ax1.twinx()

ax2.set_xlabel("Distance regressed (mm)")
ax2.set_ylabel("Burn surface (m^2)")
ax2.plot(d_reg*10**3, S_list, color='b', label='Burn Surface')
fig.tight_layout()
fig.legend()
plt.show()







