import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize
from Enginesimu import Simulation

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
alpha = 12 * np.pi / 180
eta = 4

A_t = (d_t ** 2 * np.pi / 4)


con = [d_port,d_out,l_p,alpha,eta,a,n,P_a,Pref,m_p,rho_p]
dummy = [d_port,d_out,l_p,alpha,eta,a,n,P_a,Pref,m_p,rho_p]

Index = int(input("What do you want to plot? Pressure = 1, Total impulse = 2, Mass flow = 3, Thrust = 4, Regression rate = 5 "))
IndexChange = int(input("Choose comparison: d_port = 0, d_out = 1,l_p = 2,alpha = 3,eta = 4,a = 5,n = 6,P_a = 7,Pref = 8,m_p = 9,rho_p = 10 "))
NewValue = float(input("Input new value: "))

print(con)
connew = dummy
print(con)
connew[IndexChange] = NewValue

print(con, connew)

plt.plot(Simulation(con)[0], Simulation(con)[Index], label='Old')
plt.plot(Simulation(connew)[0], Simulation(connew)[Index], label='New')
plt.xlabel("Time [s]")
plt.ylabel("Variable [SI]")
plt.grid()
plt.legend()
plt.show()

plt.title("Percentage difference")
plt.plot(Simulation(con)[0], 100 * np.abs(Simulation(con)[Index] - Simulation(con)[Index]) / Simulation(con)[Index])
plt.xlabel("Time [s]")
plt.ylabel("Difference [%]")
plt.grid()
plt.show()

# plt.plot(Simulation()[0], Simulation()[1], label='Chamber pressure')
# plt.plot(Simulation(a = 0.001)[0], Simulation(a = 0.001)[1], label='Chamber pressure eta = 2')
# plt.xlabel("Time [s]")
# plt.ylabel("Pressure [Pa]")
# plt.grid()
# plt.legend()
# plt.show()
#
# plt.plot(Simulation()[0], Simulation()[3], label='Mass flow')
# plt.xlabel("Time [s]")
# plt.ylabel("Mass flow [kg/s]")
# plt.grid()
# plt.legend()
# plt.show()
#
# plt.plot(Simulation()[0], Simulation()[4], label='Thrust')
# plt.xlabel("Time [s]")
# plt.ylabel("Thrust [N]")
# plt.grid()
# plt.legend()
# plt.show()
#
# plt.plot(Simulation()[0], Simulation()[2], label='Total impulse')
# plt.xlabel("Time [s]")
# plt.ylabel("Total impulse [Ns]")
# plt.grid()
# plt.legend()
# plt.show()

# plt.plot(Simulation()[0], Simulation()[5], label='Regression rate')
# plt.xlabel("Time [s]")
# plt.ylabel("Regression rate [m/s]")
# plt.grid()
# plt.legend()
# plt.show()
