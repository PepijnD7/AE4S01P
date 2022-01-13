import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize, linearizeAccumulate
from Read_Provided_data import plot_given_data

# TODO: Implement variation of a based on temperature variations

# Auxiliary functions
def regrate(P, a, n):
    r = a * P**n
    return r


def DivLoss(alpha):
    return (1 + np.cos(alpha)) / 2


def Kerckhove(Gamma):
    return np.sqrt(Gamma) * (2 / (1 + Gamma))**(0.5 * (Gamma + 1) / (Gamma - 1))


def FindPratio(Gamma, eps):
    Guess = 0.1
    for _ in range(10):
        Guess = (Kerckhove(Gamma)**2 / (eps**2 * (2 * Gamma / (Gamma - 1)) * (1 - Guess**((Gamma - 1) / Gamma))))**(Gamma / 2)
    return Guess


# Simulation function
def Simulation(con, density=0.0):   # Propellant density can be entered as an input, otherwise it is determined as m_p/V_grain
    g0 = 9.81
    d_port, d_out, l_p, alpha, eps, a, n, m_p, P_a, T_a = con
    a = (1.55 * 10 ** (-5) * (T_a - 273.15) + 4.47 * 10 ** (-3))* (10 ** (-6)) ** n     # Determine regression constant
                                                                                        # with ambient temp

    # INITIAL CONDITIONS and CHAMBER FILL:
    # Parameters that do not change:
    A_t = (d_t ** 2 * np.pi / 4)
    V_p = l_p * ((d_out/2)**2 - (d_port/2)**2) * np.pi

    if density == 0:
        rho_p = m_p/V_p
    else:
        rho_p = density

    # Parameters that will change:
    V_c = l_p * (d_port/2)**2 * np.pi   # Chamber volume
    P_c = P_a                           # Initial chamber pressure = ambient
    dt = 0.001                          # Small time step for chamber filling
    dpdt = 1

    m_list = []                         # Lists to append interesting data
    Isp_list = []
    p_list = []
    T_list = []
    r_list = []
    I_list = [0]
    pepa_list = []
    t_list = [0]
    m_accum = 1

    while np.abs(m_accum)>0.01:       # While the difference between m_in and m_out > 0.01, the chamber is "filling"
        S = np.pi * l_p * d_port      # Determine burning surface

        if P_c < 10**6:                 # Find the combustion data from RPA-Lite, using the right function
            Gamma = linearizeAccumulate('Gamma', P_c)
            vdk = Kerckhove(Gamma)
            c_star = linearizeAccumulate('Characteristic velocity_opt', P_c)
            C_f0 = linearizeAccumulate('Thrust coefficient_opt', P_c)
            T_c = linearizeAccumulate("Temperature", P_c)
            R = linearizeAccumulate('Gas Constant', P_c)
        else:
            Gamma = linearize('Gamma', P_c)
            vdk = Kerckhove(Gamma)
            c_star = linearize('Characteristic velocity_opt', P_c)
            C_f0 = linearize('Thrust coefficient_opt', P_c)
            T_c = linearize("Temperature", P_c)
            R = linearizeAccumulate('Gas Constant', P_c) * 1000

        rho_c = P_c/R/T_c           # Combustion gas density
        r = regrate(P_c, a, n)      # Regression rate
        m_in = r * S * rho_p        # Mass flow from grain burning
        m_out = A_t * P_c / c_star  # Mass flow through the nozzle
        m_accum = m_in - m_out      # Mass accumulating in the chamber

        dpdt = vdk**2/V_c * (c_star**2 * (rho_p - rho_c) * S * r - c_star * A_t * P_c)  # Change in chamber pressure

        C_f = C_f0 * DivLoss(alpha) + (FindPratio(Gamma, eps) - P_a / P_c) * eps    # Thrust calculation, with div. loss
        T = np.max((C_f * P_c * A_t, 0))
        Isp = T / g0 / m_out

        m_list.append(m_out)        # Adding all data to the lists
        Isp_list.append(Isp)
        p_list.append(P_c)
        pepa_list.append(P_c * FindPratio(Gamma, eps) / P_a)
        T_list.append(T)
        r_list.append(r)
        prev_I = I_list[-1]
        I_list.append(prev_I + (T * dt))
        t_list.append(t_list[-1]+dt)

        d_port+= r*dt               # Update grain port diameter and chamber pressure
        P_c += dpdt*dt
    print("CHAMBER FILLED at t=", t_list[-1],"\n")

    # OLD CHAMBER FILL
    # for _ in range(10):
    #     c_star = linearize('Characteristic velocity_opt', P_c)
    #     P_c = (c_star * rho_p * a * S / A_t)**(1 / (1-n))
    # r = regrate(P_c, a, n)
    # m_dot = r * S * rho_p
    #
    # c_star = linearize('Characteristic velocity_opt', P_c)
    # Gamma = linearize('Gamma', P_c)
    # C_f0 = linearize('Thrust coefficient_opt', P_c)
    # T_c = linearize("Temperature", P_c)
    # R = linearize('Gas constant', P_c) * 1000
    # rho_c = P_c / R / T_c
    #
    # C_f = C_f0 * DivLoss(alpha) + (FindPratio(Gamma, eps) + P_a / P_c) * eps

    # BURN LOOP
    dt=0.1                          # Use a larger time step for steady state calculation
    o = 0                           # Loop counter
    while d_port <= d_out:
        o += 1
        S = d_port * np.pi * l_p    # Determine burning surface

        P_c = (c_star * (rho_p-rho_c) * a * S / A_t) ** (1 / (1 - n))   # Update chamber pressure in steady state

        T_c = linearize("Temperature", P_c)                             # Find combustion data from RPA-Lite
        R = linearize("Gas constant", P_c) * 1000
        c_star = linearize('Characteristic velocity_opt', P_c)
        Gamma = linearize('Gamma', P_c)
        C_f0 = linearize('Thrust coefficient_opt', P_c)

        rho_c = P_c / R / T_c                                           # Update chamber parameters
        r = regrate(P_c, a, n)
        m_in = rho_p * r * S
        m_out = A_t * P_c / c_star

        C_f = C_f0 * DivLoss(alpha) + (FindPratio(Gamma, eps) - P_a / P_c) * eps    # Determine thrust
        T = C_f * P_c * A_t
        Isp = T / g0 / m_out

        m_list.append(m_out)        # Adding all data to the lists
        Isp_list.append(Isp)
        p_list.append(P_c)
        pepa_list.append(P_c*FindPratio(Gamma, eps)/P_a)
        T_list.append(T)
        r_list.append(r)
        prev_I = I_list[-1]
        I_list.append(prev_I + (T * dt))
        t_list.append(t_list[-1]+dt)

        d_port += 2 * r * dt        # Update port diameter
    print("BURN DONE AT t=", t_list[-1], "\n")

    # CHAMBER EMPTYING
    r=0                             # Regression rate = 0
    dt = 0.001
    while P_c>P_a:                  # Chamber needs to be emptied until Pc=Pa
        S = np.pi * l_p * d_port

        if P_c < 10**6:             # Get combustion data from RPA-Lite
            Gamma = linearizeAccumulate('Gamma', P_c)
            vdk = Kerckhove(Gamma)
            c_star = linearizeAccumulate('Characteristic velocity_opt', P_c)
            C_f0 = linearizeAccumulate('Thrust coefficient_opt', P_c)
            T_c = linearizeAccumulate("Temperature", P_c)
            R = linearizeAccumulate('Gas Constant', P_c)
        else:
            Gamma = linearize('Gamma', P_c)
            vdk = Kerckhove(Gamma)
            c_star = linearize('Characteristic velocity_opt', P_c)
            C_f0 = linearize('Thrust coefficient_opt', P_c)
            T_c = linearize("Temperature", P_c)
            R = linearizeAccumulate('Gas Constant', 900000)

        rho_c = P_c/R/T_c
        m_in = r * S * rho_p
        m_out = A_t * P_c / c_star
        m_accum = m_in - m_out

        dpdt = vdk**2/V_c * (c_star**2 * (rho_p - rho_c) * S * r - c_star * A_t * P_c)
        C_f = C_f0 * DivLoss(alpha) + (FindPratio(Gamma, eps) - P_a / P_c) * eps
        T = np.max((C_f * P_c * A_t, 0))
        Isp = T / g0 / m_out

        m_list.append(m_out)
        Isp_list.append(Isp)
        p_list.append(P_c)
        pepa_list.append(P_c * FindPratio(Gamma, eps) / P_a)
        T_list.append(T)
        r_list.append(r)
        prev_I = I_list[-1]
        I_list.append(prev_I + (T * dt))
        t_list.append(t_list[-1]+dt)

        P_c += dpdt*dt

    print("CHAMBER EMPTY AT t=", t_list[-1], "\n")
    I_list = I_list[1:]
    t_list = t_list[1:]

    return np.array(t_list), np.array(p_list), np.array(I_list), np.array(m_list), np.array(T_list), \
           np.array(r_list), np.array(Isp_list), np.array(pepa_list)


if __name__=='__main__':
    # set given to True for validation
    given = False
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

    const = [d_port, d_out, l_p, alphast, eps, a, n, m_p, P_a, T_a]
    conII = [d_port, d_out, l_p, alphaII, eps, a, n, m_p, P_a, T_a]
    const_given = [d_port_given, d_out_given, l_p_given, alpha_given, eps, a, n, m_p_given, P_a, T_a]
    t_st, p_st, I_st, m_st, T_st, r_st, Isp_st, pepa_st = Simulation(const)
    t_II, p_II, I_II, m_II, T_II, r_II, Isp_II, pepa_II = Simulation(conII)
    t_given, p_given, _, _, _, _, _, _ = Simulation(const_given, density=rho_p_given)

    if not given:
        # Printing output
        print("CHAMBER PRESSURE:\n")
        print("Max [MPa]:", np.max(p_II)*10**-6)
        print("Average [MPa]:", np.average(p_II)*10**-6, "\n")

        print("MASS FLOW:\n")
        print("Max:", np.max(m_II))
        print("Ave:", np.average(m_II), "\n")

        print("REGRESSION RATE:\n")
        print("Max [mm/s]:", np.max(r_II)*10**3)
        print("Ave [mm/s]:", np.average(r_II)*10**3, "\n")

        print("THRUST:\n")
        print("Max:\n 12 deg:", np.max(T_st), "\n 15 deg:", np.max(T_II), "\n Compare [%]:",
              (np.max(T_II)-np.max(T_st))/np.max(T_st)*100, "\n")

        print("Ave:\n 12 deg:", np.average(T_st), "\n 15 deg:", np.average(T_II), "\n Compare [%]:",
              (np.average(T_II)-np.average(T_st))/np.average(T_st)*100, "\n")

        print("SPECIFIC IMPULSE:\n")
        print("Max:\n 12 deg:", np.max(Isp_st), "\n 15 deg:", np.max(Isp_II), "\n Compare [%]:",
              (np.max(Isp_II)-np.max(Isp_st))/np.max(Isp_st)*100, "\n")

        print("Ave:\n 12 deg:", np.average(Isp_st), "\n 15 deg:", np.average(Isp_II), "\n Compare [%]:",
              (np.average(Isp_II)-np.average(Isp_st))/np.average(Isp_st)*100, "\n")

        print("TOTAL IMPULSE:\n")
        print("12 deg:", np.max(I_st), "\n15 deg:", np.max(I_II), "\n Compare [%]:",
              (np.max(I_II)-np.max(I_st))/np.max(I_st)*100, "\n")

        # Plotting outputs
        plt.subplot(1, 3, 1)
        plt.plot(t_st, p_II*10**-6, label='Chamber pressure')
        plt.xlabel("Time [s]", fontsize=13)
        plt.ylabel("Chamber pressure [MPa]", fontsize=13)
        plt.title("Chamber pressure", fontsize=13)
        plt.grid()

        plt.subplot(1, 3, 2)
        plt.plot(t_II, m_II, label='Mass flow')
        plt.xlabel("Time [s]", fontsize=13)
        plt.ylabel("Mass flow [kg/s]", fontsize=13)
        plt.title("Mass flow", fontsize=13)
        plt.grid()

        plt.subplot(1, 3, 3)
        plt.plot(t_II, r_II*10**3, label='Regression rate')
        plt.xlabel("Time [s]", fontsize=13)
        plt.ylabel("Regression rate [mm/s]", fontsize=13)
        plt.title("Regression rate", fontsize=13)
        plt.grid()
        plt.show()

        plt.subplot(1, 3, 1)
        plt.plot(t_II, T_II, label='Config. II')
        plt.plot(t_st, T_st, label='Standard')
        plt.xlabel("Time [s]", fontsize=13)
        plt.ylabel("Thrust [N]", fontsize=13)
        plt.title("Thrust", fontsize=13)
        plt.grid()
        plt.legend()

        plt.subplot(1, 3, 2)
        plt.plot(t_II, Isp_II, label='Config. II')
        plt.plot(t_st, Isp_st, label='Standard')
        plt.xlabel("Time [s]", fontsize=13)
        plt.ylabel("Specific impulse [s]", fontsize=13)
        plt.title("Specific impulse", fontsize=13)
        plt.grid()
        plt.legend()

        plt.subplot(1, 3, 3)
        plt.plot(t_II, I_II, label='Config. II')
        plt.plot(t_st, I_st, label='Standard')
        plt.xlabel("Time [s]", fontsize=13)
        plt.ylabel("Total impulse [Ns]", fontsize=13)
        plt.title("Total impulse", fontsize=13)
        plt.grid()
        plt.legend()
        plt.show()

        # plt.plot(t_II, pepa_II, label='Config. II')
        # plt.xlabel("Time [s]", fontsize=13)
        # plt.ylabel(r"$p_e/p_a$", fontsize=13)
        # plt.grid()
        # plt.show()


    if given:
        plt.plot(t_given, p_given * 10 ** -6, label='Chamber pressure')
        plt.xlabel("Time [s]", fontsize=13)
        plt.ylabel("Chamber pressure [MPa]", fontsize=13)
        plt.title("Chamber pressure", fontsize=13)
        plt.grid()
        plt.show()

        plot_given_data([t_given, p_given])
