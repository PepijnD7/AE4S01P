import numpy as np
from matplotlib import pyplot as plt
from Read_Data import linearize
from Read_Provided_data import plot_given_data
from read_cal_data import read_data
from read_cal_data import get_properties
#from Enginesimu import *

# Get test data of configuration from read_cal_data
Imp_test = read_data('Config2_211221_132537')['IM']
T_test = read_data('Config2_211221_132537')['LC']
Pc_test = read_data('Config2_211221_132537')['PS']
time_test = read_data('Config2_211221_132537')['LC_time']
time_Pc_test = read_data('Config2_211221_132537')['PS_time']

# Get test data of reference test
Imp_ref = read_data('ReferenceMotor_211222_092347')['IM']
T_ref = read_data('ReferenceMotor_211222_092347')['LC']
Pc_ref = read_data('ReferenceMotor_211222_092347')['PS']
time_ref = read_data('ReferenceMotor_211222_092347')['LC_time']
time_Pc_ref = read_data('ReferenceMotor_211222_092347')['PS_time']


print(get_properties('Config2_211221_132537'))
print(get_properties('ReferenceMotor_211222_092347'))

Imp_ref_list = []
T_ref_list = []
Pc_ref_list = []
Imp_test_list = []
T_test_list  = []
Pc_test_list = []
for i in range(1907,9291,1):
    Imp_test_list.append(Imp_test[i])
    T_test_list.append(T_test[i])

for i in range(5115,12499,1):
    Imp_ref_list.append(Imp_ref[i])
    T_ref_list.append(T_ref[i])

for i in range(952,4616,1):
    Pc_test_list.append(Pc_test[i])

for i in range(2556,6220,1):
    Pc_ref_list.append(Pc_ref[i])

print(len(Pc_test_list))
print(len(Pc_ref_list))
t = np.arange(0,14.768,0.002)
t_Pc = np.arange(0,14.656,0.004)

plt.subplot(1,3,1)
plt.suptitle("Comparison of reference data of original configuration with experimental data for Configuration II")
plt.plot(t,T_test_list, label = 'Experimental data config II')
plt.plot(t,T_ref_list, label = 'Reference data')
plt.ylabel("Thrust [N]")
plt.xlabel('Time [s]')
plt.grid()
plt.legend()
plt.subplot(1,3,2)
plt.plot(t_Pc,Pc_test_list,label = 'Experimental data config II')
plt.plot(t_Pc,Pc_ref_list,label = 'Reference data')
plt.ylabel('Chamber Pressure [Pa]')
plt.xlabel('Time [s]')
plt.legend()
plt.grid()
plt.subplot(1,3,3)
plt.plot(t,Imp_test_list,label = 'Experimental data config II')
plt.plot(t,Imp_ref_list,label = 'Reference data')
plt.ylabel('Total Impulse [Ns]')
plt.xlabel('Time [s]')
plt.grid()
plt.legend()
plt.show()
