import numpy as np
from matplotlib import pyplot as plt

a = 0.005132
n = 0.222
Pref = 1*10**6
Prange = np.arange(1, 500)*10**6


def regrate(P, a, n):
    r = a * (P/Pref)**n
    return r


rrange = regrate(Prange, a, n)

plt.plot(Prange, rrange)
plt.xlabel("Chamber Pressure (Pa)")
plt.ylabel("Regression rate (m/s)")
plt.grid()
plt.show()
