import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize
from Enginesimu import Simulation
from pandas import DataFrame

# Propellant data
a = 0.005132
n = 0.222

Pref = 1*10**6
P_a = 101325

m_p = 0.758
l_p = 0.107
d_out = 0.0766
d_port = 0.025

d_t = 8.37 / 1000
alpha = 12 * np.pi / 180
eps = 4
T_a = 287.15

A_t = (d_t ** 2 * np.pi / 4)

# Propellant mass variation
V_p = l_p * ((d_out/2)**2 - (d_port/2)**2) * np.pi
rho_p = m_p/V_p


con = [d_port,d_out,l_p,alpha,eps,a,n,P_a,Pref,m_p,T_a]
dummy = [d_port,d_out,l_p,alpha,eps,a,n,P_a,Pref,m_p,T_a]

Diff_P = []
Diff_I = []
Diff_mdot = []
Diff_Ft = []
Diff_r = []
Diff_Isp = []

Increment = -10
Differences = []
conSim = Simulation(con)

for i in range(len(con)):
    Difference_list_2 = []
    connew = [d_port,d_out,l_p,alpha,eps,a,n,P_a,Pref,m_p,T_a]
    connew[i] = connew[i] + Increment / 100 * connew[i]
    print(connew)

    for IO in range(1,7):
        connewSim = Simulation(connew)
        # complen = np.min((len(conSim[IO]), len(connewSim[IO])))-1
        # Difference_list = 100 * (- conSim[IO][:complen] + connewSim[IO][:complen]) / conSim[IO][:complen]
        Start_element = 100 * (connewSim[IO][0] - conSim[IO][0]) / conSim[IO][0]
        End_element = 100 * (connewSim[IO][-1] - conSim[IO][-1]) / conSim[IO][-1]
        Average_element = 100 * (np.mean(connewSim[IO]) - np.mean(conSim[IO])) / np.mean(conSim[IO])
        Difference_list_2.append([np.round(Start_element,1), np.round(End_element,1), np.round(Average_element,1)])

    Differences.append(Difference_list_2)

print("\n")
print(Differences)
print("\n")
print("d_port",Differences[0])
print("d_out",Differences[1])
print("l_p",Differences[2])
print("alpha",Differences[3])
print("eps",Differences[4])
print("a", Differences[5])
print("n",Differences[6])
print("P_a",Differences[7])
print("Pref",Differences[8])
print("Mass",Differences[9])
print("T_a",Differences[10])
print("\n")

for i in range(11):
    Diff_P.append(Differences[i][0])
    Diff_I.append(Differences[i][1])
    Diff_mdot.append(Differences[i][2])
    Diff_Ft.append(Differences[i][3])
    Diff_r.append(Differences[i][4])
    Diff_Isp.append(Differences[i][5])

df = DataFrame({'Chamber pressure': Diff_P, 'Total impulse': Diff_I, 'Mass flow': Diff_mdot , 'Thrust': Diff_Ft , 'Regression rate': Diff_r , 'Specific impulse': Diff_Isp})
print(df)

df.to_excel('SensitivityValues10D.xlsx', sheet_name='Sensitivity', index=True)

Index = int(input("What do you want to plot? Pressure = 1, Total impulse = 2, Mass flow = 3, Thrust = 4, Regression rate = 5, Specific Impulse = 6 "))
IndexChange = int(input("Choose comparison: d_port = 0, d_out = 1 , l_p = 2 , alpha = 3 , eps = 4 , a = 5 , n = 6 , P_a = 7 , Pref = 8 , m_p = 9 , T_a = 10 "))
NewValue = float(input("Percentage increase: "))

print(con)
connew = dummy
print(con)
connew[IndexChange] = connew[IndexChange] + NewValue / 100 * connew[IndexChange]

print(con, connew)

plt.plot(Simulation(con)[0], Simulation(con)[Index], label='Old')
plt.plot(Simulation(connew)[0], Simulation(connew)[Index], label='New')
plt.xlabel("Time [s]")
plt.ylabel("Variable [SI]")
plt.grid()
plt.legend()
plt.show()

plt.title("Percentage difference")
plt.plot(Simulation(con)[0], 100 * np.abs(Simulation(con)[Index] - Simulation(connew)[Index]) / Simulation(con)[Index])
plt.xlabel("Time [s]")
plt.ylabel("Difference [%]")
plt.grid()
plt.show()

# plt.plot(Simulation()[0], Simulation()[1], label='Chamber pressure')
# plt.plot(Simulation(a = 0.001)[0], Simulation(a = 0.001)[1], label='Chamber pressure eps = 2')
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