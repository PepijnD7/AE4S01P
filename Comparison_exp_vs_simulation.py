import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize
from Read_Provided_data import plot_given_data
from read_cal_data import read_data
from read_cal_data import get_properties
from Quality_factors import Simulation
from Enginesimu import dt
# from Enginesimu import Simulation

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
for i in range(0,len(time_test),int(dt/0.002)):
    time_test_list.append(time_test[i])
    T_test_list.append(T_test[i])
    Imp_test_list.append(Imp_test[i])


for i in range(0, len(Pc_test), int(dt/0.004)):
    time_Pc_list.append(time_Pc[i])
    Pc_test_list.append(Pc_test[i])


# Obtain simulation data with adapted parameters for our configuration
n = 0.222
a = 0.005132 * (10 ** (-6)) ** n
P_a = 102600
rho_p = 1720.49
m_p = 0.767
l_p = 0.10375
d_out = 0.07559
d_port = 0.02478
d_t = 0.0083
alphaII = 15 * np.pi / 180          # Configuration II
eps = 4
pepc = 0.04918
T_a = 277.15

# Quality factors
# xi_n = 1    # Nozzle quality
# xi_c = 1    # Combustion quality
# xi_d = 0.75 # Discharge coefficient

const = [d_port, d_out, d_t, l_p, alphaII, eps, a, n, m_p, P_a, T_a]
# t_II, p_II, I_II, m_II, T_II, r_II, Isp_II, pepa_II = Simulation(const, xi_n=0.855, xi_c=1.019521)
# t_II, p_II, I_II, m_II, T_II, r_II, Isp_II, pepa_II = Simulation(const, xi_n=1, xi_c=1)
t_II, p_II, I_II, m_II, T_II, r_II, Isp_II, pepa_II, gamma_list = Simulation(const, xi_n=0.88644, xi_c=1.0641235)

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
condition = abs(np.array(time_test_list) - get_properties('Config2_211221_132537', (2.125-1.390))['Start time'])
condition_pressure = abs(np.array(time_Pc_list) - get_properties('Config2_211221_132537', (2.125-1.390))['Start time'])
start_index = np.where(condition < dt)[0][0]
start_index_pressure = np.where(condition_pressure < dt)[0][0]

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
print(len(time_test_list),len(time_Pc_list))

# Compute area under pressure curves to find combustion quality
A_Pc_test = []
A_Pc_sim = []
A_t = (d_t ** 2 * np.pi / 4)
for i in range(0,len(time_Pc_list)-1):
    dA_test = A_t/ 0.735 * 0.004 * (Pc_test_list[i] + Pc_test_list[i+1])/2
    A_Pc_test.append(dA_test)
    dA_sim = A_t/0.767 * 0.004 * (Pc_sim_II_list[i] + Pc_sim_II_list[i+1])/2
    A_Pc_sim.append(dA_sim)

A_Pc_sim = sum(A_Pc_sim)
A_Pc_test = sum(A_Pc_test)
print(A_Pc_test)
print(A_Pc_sim)
print(A_Pc_test - A_Pc_sim)

print(get_properties('Config2_211221_132537', (2.125-1.390)))
# Compute thrust coefficient
gamma_avg = sum(gamma_list)/len(gamma_list)
print(gamma_avg)
T_test_list = np.array(T_test_list)
Pc_test_list = np.array(Pc_test_list)
Cf_test = 147.2/(2176649.5 * A_t) - eps * pepc + eps *  P_a/2176649.5 #np.average(T_test_list[T_test_list>2.5])/(np.average(Pc_test_list[Pc_test_list>P_a]) * A_t) - eps * pepc + eps * P_a/np.average(Pc_test_list[Pc_test_list>P_a])
Cf_sim = np.average(T_sim_II)/(np.average(Pc_sim_II) * A_t) - eps * pepc + eps *  P_a/np.average(Pc_sim_II)
print("thrust coefficient for config II is:", Cf_test)
print("thrust coefficient for simulation of config II is:",Cf_sim)
print("Combustion quality factor is:", A_Pc_test/A_Pc_sim)
print("Nozzle quality factor is:", Cf_test/Cf_sim)

# Create error plots
T_error = []
Imp_error = []
Pc_error = []
for i in range(0,len(time_test_list)):
    diff_T = abs(T_sim_II_list[i] -T_test_list[i])
    diff_Imp = abs(Imp_sim_II_list[i] - Imp_test_list[i])
    T_error.append(diff_T)
    Imp_error.append(diff_Imp)
for i in range(0,len(time_Pc_list)):
    diff_pc = abs(Pc_sim_II_list[i] - Pc_test_list[i])
    Pc_error.append(diff_pc)


line_width = 1.5
fig, ax = plt.subplots(2, 3, gridspec_kw={'height_ratios': [2, 1]})
fig.suptitle("Comparison of simulation data with test data for configuration II, \n with quality factors applied")

l1, = ax[0, 0].plot(time_test_list[500:-2200],T_test_list[500:-2200], label = 'Test data config II', linewidth=line_width)
l2, = ax[0, 0].plot(time_test_list[500:-2200],T_sim_II_list[500:-2200], label = 'Simulation data config II', linewidth=line_width)
ax[0, 0].set_ylabel("Thrust [N]")
ax[0, 0].set_xlabel('Time [s]')
ax[0,0].grid()

ax[1, 0].plot(time_test_list[500:-2200], T_error[500:-2200], linestyle="--", color='mediumaquamarine')
ax[1, 0].set_ylabel("Absolute error")
ax[1,0].grid()

ax[0, 1].plot(time_Pc_list[500:-2200],Pc_test_list[500:-2200], label='Test data config II', linewidth=line_width)
ax[0, 1].plot(time_Pc_list[500:-2200],Pc_sim_II_list[500:-2200], label='Simulation data config II', linewidth=line_width)
ax[0, 1].set_ylabel('Chamber Pressure [Pa]')
ax[0, 1].set_xlabel('Time [s]')
ax[0,1].grid()

ax[1, 1].plot(time_Pc_list[500:-2200], Pc_error[500:-2200], linestyle="--", color='mediumaquamarine')
ax[1, 1].set_ylabel("Absolute error")
ax[1,1].grid()

ax[0, 2].plot(time_test_list[500:-2200],Imp_test_list[500:-2200], label='Test data config II', linewidth=line_width)
ax[0, 2].plot(time_test_list[500:-2200],Imp_sim_II_list[500:-2200], label='Simulation data config II', linewidth=line_width)
ax[0, 2].set_ylabel('Total Impulse [Ns]')
ax[0, 2].set_xlabel('Time [s]')
ax[0,2].grid()

ax[1, 2].plot(time_test_list[500:-2200], Imp_error[500:-2200], linestyle="--", color='mediumaquamarine')
ax[1, 2].set_ylabel("Absolute error")
ax[1,2].grid()

fig.legend([l1, l2],  # The line objects
           labels=['Test data config II', 'Simulation data'],  # The labels for each line
           loc="lower center",  # Position of legend
           borderaxespad=0.1,  # Small spacing around legend box
           bbox_to_anchor=(0, 0.9, 1, 1),
           ncol=3
           )
plt.subplots_adjust(top=0.88)
plt.show()

#Printing simulation values

print("CHAMBER PRESSURE:\n")
print("Max [MPa]:", np.max(p_II) * 10 ** -6)
print("Average [MPa]:", np.average(p_II) * 10 ** -6, "\n")

print("MASS FLOW:\n")
print("Max:", np.max(m_II))
print("Ave:", np.average(m_II), "\n")

print("THRUST:\n")
print("Max:", np.max(T_II))
print("Ave: 15 deg:", np.average(T_II))

print("SPECIFIC IMPULSE:\n")
print("Ave: 15 deg:", np.average(Isp_II))

print("TOTAL IMPULSE:\n")
print("15 deg:", np.max(I_II))








