import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize
from Read_Provided_data import plot_given_data
from read_cal_data import read_data
from read_cal_data import get_properties
#from Quality_factors import Simulation
from Enginesimu import dt
from Enginesimu import Simulation

# Configuration II
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
T_a = 277.15

xi_c_II = 0.8
xi_n_II = 0.86
xi_Isp = xi_n_II * xi_c_II

# Extract simulation data for configuration II
const = [d_port, d_out, d_t, l_p, alphaII, eps, a, n, m_p, P_a, T_a]
t_II, p_II, I_II, m_II, T_II, r_II, Isp_II, pepa_II = Simulation(const)
tb_sim = t_II[-1] - t_II[0]

T_avg_sim = sum(T_II)/len(T_II)
I_tot_sim = I_II[-1]
m_avg_II_sim = sum(m_II) / len(m_II) # Average mass flow taken from list of mass flows returned from the simulation

m_avg_II_sim_b = m_p / tb_sim # Average mass flow by dividing propellant mass by burn time
Isp_avg_II_sim = T_avg_sim / (m_avg_II_sim * 9.81)  # Average Isp computed by dividing thrust by average mass flow * g0
Isp_avg_II_sim_b = I_tot_sim / (m_p * 9.81)     # Average Isp by dividing total impulse by total propellant weight

# Extract test data for configuration II
properties_II = get_properties('Config2_211221_132537')
tb = properties_II['Burn time']
I_tot = properties_II['Total Impulse']
T_max = properties_II['Max. Thrust']
T_avg = properties_II['Average Thrust']



m_avg_II = m_p / tb     #Average mass flow from test data by dividing propellant mass by burn time
Isp_avg_II = T_avg / (m_avg_II*9.81)    # Average Isp by dividing thrust by average mass flow*g0
Isp_avg_II_b = I_tot / (m_p*9.81)       # Average Isp by dividing total impulse by propellant weight

# Reference configuration

# Extract test data for reference configuration
print('Reference configuration:')

properties_ref = get_properties('ReferenceMotor_211222_092347')
tb_ref = properties_ref['Burn time']
I_tot_ref = properties_ref['Total Impulse']
T_max_ref = properties_ref['Max. Thrust']
T_avg_ref = properties_ref['Average Thrust']

m_avg_ref = m_p / tb_ref
Isp_avg_ref = T_avg_ref / (m_avg_ref*9.81)
Isp_avg_ref_b = I_tot_ref / (m_p*9.81)
print(Isp_avg_ref, Isp_avg_ref_b)


xi_c_ref = 0.83
xi_n_ref = 0.9