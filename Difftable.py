import numpy as np
from Quality_factors import Simulation

# Reference data
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
t_ref, p_ref, I_ref, m_ref, T_ref, r_II, Isp_ref, pepa_ref, gamma_list = Simulation(const, xi_n=0.92085, xi_c=0.97537)
# t_ref, p_ref, I_ref, m_ref, T_ref, r_ref, Isp_ref, pepa_ref, gamma_list = Simulation(const, xi_n=1.0, xi_c=1.0)




# Config II data
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

const = [d_port, d_out, d_t, l_p, alphaII, eps, a, n, m_p, P_a, T_a]
t_II, p_II, I_II, m_II, T_II, r_II, Isp_II, pepa_II, gamma_list = Simulation(const, xi_n=0.88644, xi_c=1.0641235)
# t_II, p_II, I_II, m_II, T_II, r_II, Isp_II, pepa_II, gamma_list = Simulation(const, xi_n=1, xi_c=1)


print("CHAMBER PRESSURE:")
print("REF Max [MPa]:", np.max(p_ref[p_ref>P_a]) * 10 ** -6)
print("REF Average [MPa]:", np.average(p_ref[p_ref>P_a]) * 10 ** -6, "\n")

print("CONFIG II Max [MPa]:", np.max(p_II[p_II>P_a]) * 10 ** -6)
print("CONFIG II Average [MPa]:", np.average(p_II[p_II>P_a]) * 10 ** -6, "\n")

print("DIFF MAX:", (np.max(p_II[p_II>P_a])-np.max(p_ref[p_ref>P_a]))/np.max(p_ref[p_ref>P_a])*100)
print("DIFF AVE:", (np.average(p_II[p_II>P_a])-np.average(p_ref[p_ref>P_a]))/np.average(p_ref[p_ref>P_a])*100, "\n")


print("THRUST:")
print("REF Max:", np.max(T_ref))
print("REF Ave:", np.average(T_ref[p_ref>P_a]), "\n")

print("CONFIG II Max:", np.max(T_II))
print("CONFIG II Ave:", np.average(T_II[p_II>P_a]), "\n")

print("DIFF MAX:", (np.max(T_II[p_II>P_a])-np.max(T_ref[p_ref>P_a]))/np.max(T_ref[p_ref>P_a])*100)
print("DIFF AVE:", (np.average(T_II[p_II>P_a])-np.average(T_ref[p_ref>P_a]))/np.average(T_ref[p_ref>P_a])*100, "\n")


print("MASS FLOW:")
print("REF Ave:", np.average(m_ref[p_ref>P_a]))
print("CONFIG II Ave:", np.average(m_II[p_II>P_a]))
print("DIFF AVE:", (np.average(m_II[p_II>P_a])- np.average(m_ref[p_ref>P_a]))/np.average(m_ref[p_ref>P_a])*100, "\n")


print("SPECIFIC IMPULSE:")
print("REF Ave:", np.average(Isp_ref[p_ref>P_a]))
print("CONFIG II Ave:", np.average(Isp_II[p_II>P_a]))
print("DIFF AVE:", (np.average(Isp_II[p_II>P_a])-np.average(Isp_ref[p_ref>P_a]))/np.average(Isp_ref[p_ref>P_a])*100, "\n")

print("TOTAL IMPULSE:")
print("REF Tot:", np.max(I_ref[p_ref>P_a]))
print("CONFIG II Tot:", np.max(I_II[p_II>P_a]))
print("DIFF Tot:", (np.max(I_II[p_II>P_a])-np.max(I_ref[p_ref>P_a]))/np.max(I_ref[p_ref>P_a])*100)
