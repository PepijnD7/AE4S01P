import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize
from Read_Provided_data import plot_given_data
from read_cal_data import read_data
from read_cal_data import get_properties
from Quality_factors import Simulation
# from Enginesimu import Simulation
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
pepc = 0.04917
T_a = 276.15

const = [d_port, d_out,d_t, l_p, alpha, eps, a, n, m_p, P_a, T_a]
t_II, p_II, I_II, m_II, T_II, r_II, Isp_II, pepa_II, gamma_list = Simulation(const,xi_n=0.92085, xi_c=0.97537)

# Create lists of simulation data
Pc_sim_II = []  # Chamber pressure during simulation (config 2)
T_sim_II = []   # Thrust during simulation (config 2)
Imp_sim_II = [] # Total impulse during simulation (config 2)

for i in range(0,len(p_II)):
    Pc_sim_II.append(p_II[i])
    T_sim_II.append(T_II[i])
    Imp_sim_II.append(I_II[i])


print(get_properties('ReferenceMotor_211222_092347', (2.158-1.424)))

# Find index of start of the burn( actual value is 40.81402 but 40.80002 is deemed accurate enough)
# Separate index is required for pressure since it has a different amount of data points
condition = abs(np.array(time_ref_list) - get_properties('ReferenceMotor_211222_092347', (2.158-1.424))['Start time'])
condition_pressure = abs(np.array(time_Pc_ref_list) - get_properties('ReferenceMotor_211222_092347', (2.158-1.424))['Start time'])
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

# Compute area under pressure curves to find combustion quality
A_Pc_ref = []
A_Pc_sim = []
A_t = (d_t ** 2 * np.pi / 4)
for i in range(0,len(time_Pc_ref_list)-1):
    dA_test = A_t/ 0.734 * 0.004 * (Pc_ref_list[i] + Pc_ref_list[i+1])/2
    A_Pc_ref.append(dA_test)
    dA_sim = A_t/0.775 * 0.004 * (Pc_sim_II_list[i] + Pc_sim_II_list[i+1])/2
    A_Pc_sim.append(dA_sim)

A_Pc_sim = sum(A_Pc_sim)
A_Pc_ref = sum(A_Pc_ref)
print(A_Pc_ref)
print(A_Pc_sim)
print(A_Pc_ref - A_Pc_sim)

# Compute thrust coefficient
gamma_avg  = sum(gamma_list)/len(gamma_list)
print(gamma_avg)
T_ref_list = np.array(T_ref_list)
Pc_ref_list = np.array(Pc_ref_list)
Cf_ref = 148.7/(2114340.7 * A_t) - eps * pepc + eps *  P_a/2114340.7 #np.average(T_ref_list[T_ref_list>2.5])/(np.average(Pc_ref_list[Pc_ref_list>P_a]) * A_t) - eps * pepc + eps * P_a/np.average(Pc_ref_list[Pc_ref_list>P_a])
Cf_sim = np.average(T_sim_II)/(np.average(Pc_sim_II) * A_t) - eps * pepc + eps *  P_a/np.average(Pc_sim_II)
print("thrust coefficient for reference config is:", Cf_ref)
print("thrust coefficient for simulation of reference config is:",Cf_sim)
print("Combustion quality factor is:", A_Pc_ref/A_Pc_sim)
print("Nozzle quality factor is:", Cf_ref/Cf_sim)


line_width = 1.5
fig, ax = plt.subplots(2, 3, gridspec_kw={'height_ratios': [2, 1]})
fig.suptitle("Comparison of simulation data with reference data from day 2 \n with quality factors applied")

l1, = ax[0, 0].plot(time_ref_list[2000:-500],T_ref_list[2000:-500], label = 'Reference data', linewidth=line_width)
l2, = ax[0, 0].plot(time_ref_list[2000:-500],T_sim_II_list[2000:-500], label = 'Simulation data reference config', linewidth=line_width)
ax[0, 0].set_ylabel("Thrust [N]")
ax[0, 0].set_xlabel('Time [s]')
ax[0,0].grid()

ax[1, 0].plot(time_ref_list[2000:-500], T_error[2000:-500], linestyle="--", color='mediumaquamarine')
ax[1, 0].set_ylabel("Absolute error [N]")
ax[1,0].grid()


ax[0, 1].plot(time_Pc_ref_list[2000:-500],Pc_ref_list[2000:-500],label = 'Reference data', linewidth=line_width)
ax[0, 1].plot(time_Pc_ref_list[2000:-500],Pc_sim_II_list[2000:-500],label = 'Simulation data reference config', linewidth=line_width)
ax[0, 1].set_ylabel('Chamber Pressure [Pa]')
ax[0, 1].set_xlabel('Time [s]')
ax[0,1].grid()

ax[1, 1].plot(time_Pc_ref_list[2000:-500], Pc_error[2000:-500], linestyle="--", color='mediumaquamarine')
ax[1, 1].set_ylabel("Absolute error [Pa]")
ax[1,1].grid()

ax[0, 2].plot(time_ref_list[2000:-500],Imp_ref_list[2000:-500],label = 'Reference data', linewidth=line_width)
ax[0, 2].plot(time_ref_list[2000:-500],Imp_sim_II_list[2000:-500],label = 'Simulation data reference config', linewidth=line_width)
ax[0, 2].set_ylabel('Total Impulse [Ns]')
ax[0, 2].set_xlabel('Time [s]')
ax[0,2].grid()

ax[1, 2].plot(time_ref_list[2000:-500], Imp_error[2000:-500], linestyle="--", color='mediumaquamarine')
ax[1, 2].set_ylabel("Absolute error [Ns]")
ax[1,2].grid()

fig.legend([l1, l2],  # The line objects
           labels=['Test data config II', 'Simulation data'],  # The labels for each line
           loc="lower center",  # Position of legend
           borderaxespad=0.1,  # Small spacing around legend box
           bbox_to_anchor=(0, 0.9, 1, 1),
           ncol=2
           )
plt.subplots_adjust(top=0.88)
plt.show()


#Printing simulation values

print("CHAMBER PRESSURE:")
print("Max [MPa]:", np.max(p_II) * 10 ** -6)
print("Average [MPa]:", np.average(p_II) * 10 ** -6, "\n")

print("MASS FLOW:")
print("Max:", np.max(m_II))
print("Ave:", np.average(m_II), "\n")

print("THRUST:")
print("Max:", np.max(T_II))
print("Ave: 12 deg:", np.average(T_II), "\n")

print("SPECIFIC IMPULSE:")
print("Ave: 12 deg:", np.average(Isp_II), "\n")

print("TOTAL IMPULSE:")
print("12 deg:", np.max(I_II))

