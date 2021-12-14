import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize

# Propellant data
a = 0.005132
n = 0.222
Pref = 1*10**6
m_p = 0.758
l_p = 0.107
d_out = 0.0766
d_port = 0.025
d_t = 8.37 / 1000
alpha = 15 * np.pi / 180


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
    d_list = d_ilist + 2*r*t
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
# plt.show()

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


# Initial conditions
P_c = 5 * Pref
S = np.pi * l_p * d_port
A_t = (d_t ** 2 * np.pi / 4)
a *= (10 ** (-6)) ** n

for _ in range(10):
    c_star = linearize('Characteristic velocity_opt', P_c)
    # print('c_star: ', c_star)

    P_c = (c_star * rho_P * a * S / A_t)**(1 / (1-n))
    # print('P_c: ', P_c)

a = 0.005132
r = regrate(P_c, a, n)
m_dot = r * S * rho_P

c_star = linearize('Characteristic velocity_opt', P_c)
C_f = linearize('Thrust coefficient_opt', P_c)
lamda = (1 + np.cos(alpha)) / 2
T = C_f * P_c * A_t * lamda

print(P_c, r, m_dot, T)


# Simulation
dt = 0.1

m_list = []
p_list = []
T_list = []
r_list = []
I_list = [0]

while d_port <= d_out:
    S_burn = d_port * np.pi * l_p

    a *= (10 ** (-6)) ** n
    P_c = (c_star * rho_P * a * S_burn / A_t) ** (1 / (1 - n))

    a = 0.005132
    r_rate = regrate(P_c, a, n)

    m_dot = rho_P * r_rate * S_burn
    c_star = linearize('Characteristic velocity_opt', P_c)
    T = m_dot * c_star * linearize('Thrust coefficient_opt', P_c)

    m_list.append(m_dot)
    p_list.append(P_c)
    T_list.append(T)
    r_list.append(r_rate)
    prev_I = I_list[-1]
    I_list.append(prev_I + (T * dt))

    d_port += 2 * r_rate * dt


t_list = np.arange(0, dt*len(p_list), dt)
I_list = I_list[1:]

plt.plot(t_list, p_list, label='Chamber pressure')
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.grid()
plt.legend()
plt.show()

plt.plot(t_list, m_list, label='Mass flow')
plt.xlabel("Time [s]")
plt.ylabel("Mass flow [kg/s]")
plt.grid()
plt.legend()
plt.show()

plt.plot(t_list, T_list, label='Thrust')
plt.xlabel("Time [s]")
plt.ylabel("Thrust [N]")
plt.grid()
plt.legend()
plt.show()

plt.plot(t_list, I_list, label='Total impulse')
plt.xlabel("Time [s]")
plt.ylabel("Total impulse [Ns]")
plt.grid()
plt.legend()
plt.show()

plt.plot(t_list, r_list, label='Regression rate')
plt.xlabel("Time [s]")
plt.ylabel("Regression rate [m/s]")
plt.grid()
plt.legend()
plt.show()
