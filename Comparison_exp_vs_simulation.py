import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize
from Read_Provided_data import plot_given_data
from read_cal_data import read_data
from read_cal_data import get_properties
from Enginesimu import *

# Get data from read_cal_data
Imp_test = read_data('Config2_211221_132537')['IM']
T_test = read_data('Config2_211221_132537')['LC']
Pc_test = read_data('Config2_211221_132537')['PS']
time_test = read_data('Config2_211221_132537')['LC_time']
time_Pc = read_data('Config2_211221_132537')['PS_time']

# Make empty lists in which data will be put with time interval of 0.1s, corresponding with simulation data
time_test_list = [] # Timestamps of the test data
time_Pc_list = []   # Timestamps of the pressure test data
Pc_test_list = []   # Pressure measurements of test data
Imp_test_list = []  # Total impulse measurements of test data
T_test_list = []    # Thrust measurements of test data

# Create lists of measurements with a time interval of 0.1 seconds
for i in range(0,len(time_test),50):
    time_test_list.append(time_test[i])
    T_test_list.append(T_test[i])
    Imp_test_list.append(Imp_test[i])

for i in range(0,len(Pc_test),25):
    time_Pc_list.append(time_Pc[i])
    Pc_test_list.append(Pc_test[i])

# Create lists of simulation data
Pc_sim_II = []  # Chamber pressure during simulation (config 2)
T_sim_II = []   # Thrust during simulation (config 2)
Imp_sim_II = [] # Total impulse during simulation (config 2)

for i in range(0,len(p_II)):
    Pc_sim_II.append(p_II[i])
    T_sim_II.append(T_II[i])
    Imp_sim_II.append(I_II[i])

# Find index of start of the burn( actual value is 40.81402 but 40.80002 is deemed accurate enough)
# Separate index is required for pressure since it has a different amount of data points
start_index = np.where(np.array(time_test_list)==40.80002)[0][0]
start_index_pressure = np.where(np.array(time_Pc_list) == 40.80004)[0][0]

T_sim_II_list = [0]*start_index     # List of zeros prior to simulation data
Imp_sim_II_list = [0]*start_index   # List of zeros prior to simulation data
zero_after = [0]*(len(time_test_list) - start_index - len(T_sim_II))                    # List of zeros after simulation
Impulse_after = [Imp_sim_II[-1]]*(len(time_test_list) - start_index - len(Imp_sim_II))  # List of Final total impulse after simulation

# Create list of Thrust and Impulse data compatible with test data with dt = 0.1s
for i in range(0,len(p_II)):
    T_sim_II_list.append(T_sim_II[i])
    Imp_sim_II_list.append((Imp_sim_II[i]))
for i in range(0,(len(time_test_list) - len(T_sim_II_list))):
    T_sim_II_list.append(zero_after[i])
    Imp_sim_II_list.append((Impulse_after[i]))

Pc_sim_II_list = [0]*start_index_pressure                                       # List of zeros prior to simulation
zero_after = [0]*(len(time_Pc_list) - start_index_pressure - len(Pc_sim_II))    # List of zeros after simulation

# Create list of Pressure data compatible with test data with dt = 0.1 s
# !!! Note, these are 249 entries instead of 250
for i in range(0,len(p_II)):
    Pc_sim_II_list.append(Pc_sim_II[i])
for i in range(0,(len(time_Pc_list) - len(Pc_sim_II_list))):
    Pc_sim_II_list.append(zero_after[i])

plt.subplot(1,3,1)
plt.suptitle("Comparison of simulation data with experimental data for Configuration II")
plt.plot(time_test_list,T_test_list, label = 'Experimental data config II')
plt.plot(time_test_list,T_sim_II_list, label = 'Simulation data config II')
plt.ylabel("Thrust [N]")
plt.xlabel('Time [s]')
plt.grid()
plt.legend()
plt.subplot(1,3,2)
plt.plot(time_Pc_list,Pc_test_list,label = 'Experimental data config II')
plt.plot(time_Pc_list,Pc_sim_II_list,label = 'Simulation data config II')
plt.ylabel('Chamber Pressure [Pa]')
plt.xlabel('Time [s]')
plt.legend()
plt.grid()
plt.subplot(1,3,3)
plt.plot(time_test_list,Imp_test_list,label = 'Experimental data config II')
plt.plot(time_test_list,Imp_sim_II_list,label = 'Simulation data config II')
plt.ylabel('Total Impulse [Ns]')
plt.xlabel('Time [s]')
plt.grid()
plt.legend()
plt.show()


print(get_properties('Config2_211221_132537'))




