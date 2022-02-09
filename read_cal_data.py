#from nptdms import TdmsFile
import numpy as np
import matplotlib.pyplot as plt


def read_tdms(filename, param_list):
    data_dict = {}

    with TdmsFile.open(filename) as tdms_file:
        for par in param_list:
            try:
                data_object = tdms_file[par + ' data'][par + '0']
            except KeyError:
                data_object = tdms_file[par + ' data'][par + '1']

            t_step  = float(data_object.properties['wf_increment'])
            t_start = float(data_object.properties['wf_start_offset'])
            t_end   = float(data_object.properties['wf_samples']) * t_step + t_start - t_step/2
            t_list  = np.arange(t_start, t_end, t_step)

            data_dict[par] = [t_list, np.array(data_object)]

    return data_dict


def calibration():
    # Read calibration file
    param_list = ['LC', 'PS']
    cal_original = read_tdms('test_data_2021/Calibration_211217_105057.tdms', param_list)

    # 'filter' values
    res = 1000
    delta = 20
    max_diffs = [0.3 * 10**(-5), 0.3 * 10**(-4)]
    ranges_list = [[[-0.00007, 0.00001], [-0.0001, 0.00001], [-0.00021, 0.00002]], [[0.004, 0.0001], [0.0049, 0.0001]]]

    # 'calibration' values
    values = [[0, 3.938 * 9.80665, 17.593 * 9.80665], [0, 7 * 10**5]]  # N, Pa
    slope_list = []
    offset_list = []

    for param, max_diff, ranges, val in zip(param_list, max_diffs, ranges_list, values):
        # 'filter' the data
        # cal_t  = cal_original[param][0][::res]
        cal_data = cal_original[param][1][::res]

        diff = abs(cal_data[:-delta] - cal_data[delta:])
        filtered = cal_data[np.where(diff < max_diff)]

        # find the 'equilibrium point' for each calibration value
        avg_list = [np.average(filtered[np.where(abs(filtered - r[0]) < r[1])]) for r in ranges]

        # plot the equilibrium points as horizontal lines over the data
        # plt.plot(cal_t, cal_data)
        # # plt.scatter(cal_t[np.where(diff < max_diff)], filtered, s=3)
        # for avg in avg_list:
        #     plt.hlines(avg, cal_t[0], cal_t[-1], colors='r')
        # plt.show()

        # find the 'slopes' and 'offsets' of the calibration
        slope1 = (val[1]  - val[0])  / (avg_list[1]  - avg_list[0])
        slope2 = (val[-1] - val[0])  / (avg_list[-1] - avg_list[0])
        slope3 = (val[-1] - val[-2]) / (avg_list[-1] - avg_list[-2])
        slope_list.append((slope1 + slope2 + slope3) / 3)

        offset_list.append(avg_list[0])

    return slope_list, offset_list


def convert(filename, begin):
    params = ['LC', 'PS']
    filtered_dict = {}

    # get calibration values and read in data
    slopes, offsets = calibration()
    data = read_tdms(filename, params)

    for param, s, off in zip(params, slopes, offsets):
        # print(data[param][1][0], off)

        # slice data to time interval where motor is fired, reduce #data-points by factor 100, subtract offset from data
        t_range = np.where(np.logical_and(data[param][0] > begin, data[param][0] < begin + 25))
        data_t = data[param][0][t_range][::100].round(5)
        data_y = data[param][1][t_range][::100] - off

        # filter data by averaging over 'n_points' number of points
        n_points = 30
        if param == 'LC':
            n_points = 1
            data_y += off - data[param][1][0]  # dit is gefoefel voor reference fire van dag 2, maar Jeije approved

        # filter the data
        filtered = data_y.copy()
        for i in range(n_points):
            filtered[:-n_points] += data_y[i:-n_points + i]
        filtered /= (n_points + 1)

        # make impulse data
        if param == 'LC':
            dt = data_t[1] - data_t[0]
            impulse = np.zeros(len(filtered[:-n_points]) + 1)

            for j in range(1, len(impulse)):
                impulse[j] = impulse[j-1] + filtered[j-1] * dt
            impulse = impulse[1:]

            # multiply data by slope from calibration and plot
            filtered_dict['IM'] = [data_t[:-n_points], (impulse * s).round(3)]
            plt.plot(*filtered_dict['IM'])
            plt.show()

        # multiply data by slope from calibration and plot
        filtered_dict[param] = [data_t[:-n_points], (filtered[:-n_points] * s).round(3)]
        plt.plot(*filtered_dict[param])
        plt.show()

    return filtered_dict


def save_data(filename_and_begin_list):
    with open('Filtered_Data.txt', 'w') as file:
        for [filename, t_begin] in filename_and_begin_list:
            file.write(filename + '\n')
            filtered_data = convert(filename, t_begin)

            for pa in filtered_data.keys():
                file.write(pa + '_time: ' + str(filtered_data[pa][0].tolist())[1:-1] + '\n')
                file.write(pa + ': '      + str(filtered_data[pa][1].tolist())[1:-1] + '\n')


def read_data(filename, corrected=True):
    with open('Filtered_Data.txt') as file:
        first_line = np.infty
        data_dict = {}

        for i, line in enumerate(file):
            if filename in line.strip():
                first_line = i + 1

            if i >= first_line + 6:
                break

            if i >= first_line:
                param, data = line.strip().split(': ')
                data_dict[param] = np.array(data.split(', '), dtype=float)

    if corrected:
        # fig, (axis1, axis2) = plt.subplots(1, 2)
        # axis1.plot(data_dict['LC_time'], data_dict['LC'], label='not corrected', c='blue')
        # axis2.plot(data_dict['IM_time'], data_dict['IM'], label='not corrected', c='blue')

        # get properties
        max_T = np.max(data_dict['LC'])
        max_T_ind = np.where(data_dict['LC'] == max_T)[0][0]
        end_T = np.average(data_dict['LC'][-100:])

        # get correction factor
        delta_Thrust = max_T - end_T
        factor_corr = max_T / delta_Thrust

        # correct thrust data
        data_dict['LC'][max_T_ind:] -= end_T
        data_dict['LC'][max_T_ind:] *= factor_corr

        # recalculate impulse
        dt = data_dict['LC_time'][1] - data_dict['LC_time'][0]
        impulse = np.zeros(len(data_dict['LC']) + 1)

        for j in range(1, len(impulse)):
            impulse[j] = impulse[j - 1] + data_dict['LC'][j - 1] * dt
        data_dict['IM'] = impulse[1:]

        # axis1.plot(data_dict['LC_time'][max_T_ind:], data_dict['LC'][max_T_ind:], label='corrected', c='r')
        # axis2.plot(data_dict['IM_time'][max_T_ind:], data_dict['IM'][max_T_ind:], label='corrected', c='r')
        # plt.show()

    return data_dict


def get_properties(filename, m_prop, corrected=True):
    data = read_data(filename, corrected)
    par_list = [['IM_time', 'IM'], ['LC_time', 'LC'], ['PS_time', 'PS']]
    par_dict = {'IM': 'Impulse', 'LC': 'Thrust', 'PS': 'Chamber Pressure'}
    prop_dict = {}

    # find begin and end of the burn
    start_value = np.average(data['LC'][:100]) + 5  # 5 is added due to noise
    ignit_index = np.where(data['LC'] > 10)[0][0]
    ignition = data['LC_time'][ignit_index]
    start = data['LC_time'][np.where(data['LC'][ignit_index + 500:] > start_value)[0][0] + ignit_index + 350]

    max_index = np.where(data['LC'] == np.max(data['LC']))[0][0]
    prop_dict['max_LC_index'] = max_index

    end_value = np.average(data['LC'][-100:]) + 5  # 5 is added due to noise
    end_index = np.where(data['LC'][max_index:] < end_value)[0][0] + max_index
    end = data['LC_time'][end_index]

    prop_dict['Burn time']   = (end - start).round(3)
    prop_dict['Ignition time'] = (start - ignition).round(3)
    prop_dict['Start time'] = start.round(3)
    prop_dict['LC_end_value'] = (end_value - 5).round(5)  # 5 is added due to noise

    # analyse data
    for [t, y] in par_list:
        # plot data with begin and end indicated as vertical line
        # plt.plot(data[t], data[y])
        # plt.vlines([ignition, start, end], np.min(data[y]), np.max(data[y]), colors='r')
        # plt.show()

        # firing properties
        if y == 'IM':
            prop_dict['Total ' + par_dict[y]] = data[y][np.where(data[t] - end == np.min(np.abs(data[t] - end)))[0][0]].round(1)

        else:
            closest_start = np.where((data[t] - start).round(4) == np.min(np.abs(data[t] - start)).round(4))[0][0]
            closest_end   = np.where((data[t] - end).round(4)   == np.min(np.abs(data[t] - end)).round(4))[0][0]

            prop_dict['Max. ' + par_dict[y]] = np.max(data[y]).round(1)
            prop_dict['Avg. ' + par_dict[y]] = np.average(data[y][closest_start:closest_end]).round(1)

    prop_dict['Avg. Mass flow'] = (m_prop / prop_dict['Burn time']).round(3)
    prop_dict['Avg. Isp'] = (prop_dict['Total Impulse'] / (9.80665 * m_prop)).round(3)

    return prop_dict


if __name__ == '__main__':
    # convert('test_data_2021/Config2_211221_132537.tdms', 34)
    # convert('test_data_2021/ReferenceMotor_211221_114135.tdms', 222)
    # convert('test_data_2021/ReferenceMotor_211222_092347.tdms', 65)

    # save_data([['test_data_2021/Config2_211221_132537.tdms', 34],
    #            ['test_data_2021/ReferenceMotor_211221_114135.tdms', 222],
    #            ['test_data_2021/ReferenceMotor_211222_092347.tdms', 65]])

    # data1 = read_data('Config2_211221_132537')
    # data2 = read_data('ReferenceMotor_211221_114135')
    # data3 = read_data('ReferenceMotor_211222_092347')

    # set 'True' to 'False' for data without strain correction
    print(get_properties('Config2_211221_132537', (2.125-1.390), True))
    print(get_properties('ReferenceMotor_211221_114135', (2.158-1.424), True))
    print(get_properties('ReferenceMotor_211222_092347', (2.158-1.424), True))  # same as ref from day 1, since we don't have it
