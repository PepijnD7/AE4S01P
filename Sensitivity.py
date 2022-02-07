import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize
from Enginesimu2 import Simulation
from pandas import DataFrame

# Propellant data
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
alphast = 12 * np.pi / 180  # Standard setup
alphaII = 15 * np.pi / 180  # Configuration II
eps = 4
T_a = 273.15 + 15

# Given data
l_p_given = 99.60 / 1000
d_out_given = 76.6 / 1000
d_port_given = 24.90 / 1000
m_p_given = 741.8 / 1000
alpha_given = 12 * np.pi / 180
rho_p_given = 1786.5


con = [d_port, d_out, d_t, l_p, alphast, eps, a, n, m_p, P_a, T_a]
dummy = [d_port, d_out, d_t, l_p, alphast, eps, a, n, m_p, P_a, T_a]

Diff_P = []
Diff_I = []
Diff_mdot = []
Diff_Ft = []
Diff_r = []
Diff_Isp = []

Increment = -10
Differences = []
conSim = Simulation(con)

# d_port, d_out, d_t, l_p, alpha, eps, a, n, m_p, P_a, T_a

for i in range(len(con)):
    Difference_list_2 = []
    connew = [d_port, d_out, d_t, l_p, alphast, eps, a, n, m_p, P_a, T_a]
    connew[i] = connew[i] + Increment / 100 * connew[i]
    print(connew)

    for IO in range(0,8):
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
print("d_t",Differences[2])
print("l_p",Differences[3])
print("alpha",Differences[4])
print("eps",Differences[5])
print("a", Differences[6])
print("n",Differences[7])
print("Mass",Differences[8])
print("P_a",Differences[9])
print("T_a",Differences[10])
print("\n")

for i in range(11):
    Diff_P.append(Differences[i][1])
    Diff_I.append(Differences[i][2])
    Diff_mdot.append(Differences[i][3])
    Diff_Ft.append(Differences[i][4])
    Diff_r.append(Differences[i][5])
    Diff_Isp.append(Differences[i][6])

df = DataFrame({'Chamber pressure': Diff_P, 'Total impulse': Diff_I, 'Mass flow': Diff_mdot , 'Thrust': Diff_Ft , 'Regression rate': Diff_r , 'Specific impulse': Diff_Isp})
print(df)

df.to_excel('SensitivityValues10D.xlsx', sheet_name='Sensitivity', index=True)

