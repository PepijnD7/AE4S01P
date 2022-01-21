import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize
from Read_Provided_data import plot_given_data
from read_cal_data import read_data
from read_cal_data import get_properties
from Quality_factors import Simulation
#from Enginesimu import Simulation
from Enginesimu import dt

# Get test data of reference test
Imp_ref = read_data('ReferenceMotor_211222_092347')['IM']
T_ref = read_data('ReferenceMotor_211222_092347')['LC']
Pc_ref = read_data('ReferenceMotor_211222_092347')['PS']
time_ref = read_data('ReferenceMotor_211222_092347')['LC_time']
time_Pc_ref = read_data('ReferenceMotor_211222_092347')['PS_time']

Imp_ref_list = []
T_ref_list = []
Pc_ref_list = []
time_ref_list = []
time_Pc_ref_list = []

# Create lists of measurements with a time interval of 0.1 seconds
for i in range(0,len(time_ref),int(dt/0.002)):
    time_ref_list.append(time_ref[i])
    T_ref_list.append(T_ref[i])
    Imp_ref_list.append(Imp_ref[i])

if dt > 0.004:
    for i in range(0,len(Pc_ref),int(dt/0.004)):
        time_Pc_ref_list.append(time_Pc_ref[i])
        Pc_ref_list.append(Pc_ref[i])
else:
    for i in range(0,len(Pc_ref),1):
        time_Pc_ref_list.append(time_Pc_ref[i])
        Pc_ref_list.append(Pc_ref[i])

# Obtain simulation data with adapted parameters for reference test configuration
n = 0.222
a = 0.005132 * (10 ** (-6)) ** n
P_a = 102965
rho_p = 1720.49
m_p = 0.775
l_p = 0.1041
d_out = 0.0755
d_port = 0.0236
d_t = 8.3 / 1000
alpha = 12 * np.pi / 180          # reference Configuration
eps = 4
T_a = 276.15

const = [d_port, d_out,d_t, l_p, alpha, eps, a, n, m_p, P_a, T_a]
t_II, p_II, I_II, m_II, T_II, r_II, Isp_II, pepa_II = Simulation(const)

# Create lists of simulation data
Pc_sim_II = []  # Chamber pressure during simulation (config 2)
T_sim_II = []   # Thrust during simulation (config 2)
Imp_sim_II = [] # Total impulse during simulation (config 2)

for i in range(0,len(p_II)):
    Pc_sim_II.append(p_II[i])
    T_sim_II.append(T_II[i])
    Imp_sim_II.append(I_II[i])


print(get_properties('ReferenceMotor_211222_092347'))

# Find index of start of the burn( actual value is 40.81402 but 40.80002 is deemed accurate enough)
# Separate index is required for pressure since it has a different amount of data points
condition = abs(np.array(time_ref_list) - get_properties('ReferenceMotor_211222_092347')['Start time'])
condition_pressure = abs(np.array(time_Pc_ref_list) - get_properties('ReferenceMotor_211222_092347')['Start time'])
start_index = np.where(condition < dt)[0][0]
start_index_pressure = np.where(condition_pressure < dt)[0][0]

T_sim_II_list = [0]*start_index     # List of zeros prior to simulation data
Imp_sim_II_list = [0]*start_index   # List of zeros prior to simulation data
zero_after = [0]*(len(time_ref_list) - start_index - len(T_sim_II))                    # List of zeros after simulation
Impulse_after = [Imp_sim_II[-1]]*(len(time_ref_list) - start_index - len(Imp_sim_II))  # List of Final total impulse after simulation

# Create list of Thrust and Impulse data compatible with test data with dt = 0.1s
for i in range(0,len(p_II)):
    T_sim_II_list.append(T_sim_II[i])
    Imp_sim_II_list.append((Imp_sim_II[i]))
for i in range(0,(len(time_ref_list) - len(T_sim_II_list))):
    T_sim_II_list.append(zero_after[i])
    Imp_sim_II_list.append((Impulse_after[i]))

Pc_sim_II_list = [0]*start_index_pressure                                       # List of zeros prior to simulation
zero_after = [0]*(len(time_Pc_ref_list) - start_index_pressure - len(Pc_sim_II))    # List of zeros after simulation

# Create list of Pressure data compatible with test data with dt = 0.1 s
# !!! Note, these are 249 entries instead of 250
for i in range(0,len(p_II)):
    Pc_sim_II_list.append(Pc_sim_II[i])
for i in range(0,(len(time_Pc_ref_list) - len(Pc_sim_II_list))):
    Pc_sim_II_list.append(zero_after[i])

T_error = []
Imp_error = []
Pc_error = []
for i in range(0,len(time_ref_list)):
    diff_T = abs(T_sim_II_list[i] -T_ref_list[i])
    diff_Imp = abs(Imp_sim_II_list[i] - Imp_ref_list[i])
    T_error.append(diff_T)
    Imp_error.append(diff_Imp)
for i in range(0,len(time_Pc_ref_list)):
    diff_pc = abs(Pc_sim_II_list[i] - Pc_ref_list[i])
    Pc_error.append(diff_pc)


plt.subplot(1,3,1)
plt.suptitle("Comparison of simulation data with reference data from day 2 \n without quality factors")
plt.plot(time_ref_list[2000:-500],T_ref_list[2000:-500], label = 'Reference data')
plt.plot(time_ref_list[2000:-500],T_sim_II_list[2000:-500], label = 'Simulation data reference config')
plt.ylabel("Thrust [N]")
plt.xlabel('Time [s]')
plt.grid()
plt.legend()
plt.subplot(1,3,2)
plt.plot(time_Pc_ref_list[2000:-500],Pc_ref_list[2000:-500],label = 'Reference data')
plt.plot(time_Pc_ref_list[2000:-500],Pc_sim_II_list[2000:-500],label = 'Simulation data reference config')
plt.ylabel('Chamber Pressure [Pa]')
plt.xlabel('Time [s]')
plt.legend()
plt.grid()
plt.subplot(1,3,3)
plt.plot(time_ref_list[2000:-500],Imp_ref_list[2000:-500],label = 'Reference data')
plt.plot(time_ref_list[2000:-500],Imp_sim_II_list[2000:-500],label = 'Simulation data reference config')
plt.ylabel('Total Impulse [Ns]')
plt.xlabel('Time [s]')
plt.grid()
plt.legend()
plt.show()

plt.subplot(1,3,1)
plt.suptitle("Difference between test data and simulation data of the reference configuration")
plt.plot(time_ref_list[2000:-500],T_error[2000:-500], label = 'Difference in thrust')
plt.ylabel("Thrust [N]")
plt.xlabel('Time [s]')
plt.grid()
plt.legend()
plt.subplot(1,3,2)
plt.plot(time_Pc_ref_list[2000:-500],Pc_error[2000:-500],label = 'Difference in chamber pressure')
plt.ylabel('Chamber Pressure [Pa]')
plt.xlabel('Time [s]')
plt.legend()
plt.grid()
plt.subplot(1,3,3)
plt.plot(time_ref_list[2000:-500],Imp_error[2000:-500],label = 'Difference in total impulse')
plt.ylabel('Total Impulse [Ns]')
plt.xlabel('Time [s]')
plt.grid()
plt.legend()
plt.show()