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

T_error = []
Imp_error = []
Pc_error = []
for i in range(0,len(Imp_ref_list)):
    diff_T = abs(T_test_list[i] -T_ref_list[i])
    diff_Imp = abs(Imp_test_list[i] - Imp_ref_list[i])
    T_error.append(diff_T)
    Imp_error.append(diff_Imp)
for i in range(0,len(Pc_test_list)):
    diff_pc = abs(Pc_test_list[i] - Pc_ref_list[i])
    Pc_error.append(diff_pc)

line_width = 1.5
fig, ax = plt.subplots(2, 3, gridspec_kw={'height_ratios': [2, 1]})
fig.suptitle("Comparison of reference data of day 2 with test data for Configuration II")

l1, = ax[0, 0].plot(t[1000:-1000],T_test_list[1000:-1000], label = 'Test data config II', linewidth=line_width)
l2, = ax[0, 0].plot(t[1000:-1000],T_ref_list[1000:-1000], label = 'Reference data', linewidth=line_width)
ax[0, 0].set_ylabel("Thrust [N]")
ax[0, 0].set_xlabel('Time [s]')

ax[1, 0].plot(t[1000:-1000], T_error[1000:-1000], linestyle="--", color='mediumaquamarine')
ax[1, 0].set_ylabel("Absolute error [N]")

ax[0, 1].plot(t_Pc[:-1000],Pc_test_list[:-1000],label = 'Test data config II', linewidth=line_width)
ax[0, 1].plot(t_Pc[:-1000],Pc_ref_list[:-1000],label = 'Reference data', linewidth=line_width)
ax[0, 1].set_ylabel('Chamber Pressure [Pa]')
ax[0, 1].set_xlabel('Time [s]')

ax[1, 1].plot(t_Pc[:-1000], Pc_error[:-1000], linestyle="--", color='mediumaquamarine')
ax[1, 1].set_ylabel("Absolute error [Pa]")

ax[0, 2].plot(t[1000:-1000],Imp_test_list[1000:-1000],label = 'Test data config II', linewidth=line_width)
ax[0, 2].plot(t[1000:-1000],Imp_ref_list[1000:-1000],label = 'Reference data', linewidth=line_width)
ax[0, 2].set_ylabel('Total Impulse [Ns]')
ax[0, 2].set_xlabel('Time [s]')

ax[1, 2].plot(t[1000:-1000], Imp_error[1000:-1000], linestyle="--", color='mediumaquamarine')
ax[1, 2].set_ylabel("Absolute error [Ns]")

fig.legend([l1, l2],  # The line objects
           labels=['Test data config II', 'Reference data'],  # The labels for each line
           loc="lower center",  # Position of legend
           borderaxespad=0.1,  # Small spacing around legend box
           bbox_to_anchor=(0, 0.83, 1, 1),
           ncol=2
           )
plt.subplots_adjust(top=0.80)
plt.show()


