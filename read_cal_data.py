from nptdms import TdmsFile
import numpy as np
import matplotlib.pyplot as plt


def read_tdms(filename, param_list):
    data_dict = {}

    with TdmsFile.open(filename) as tdms_file:
        for par in param_list:
            data_object = tdms_file[par + ' data'][par + '0']

            t_step  = float(data_object.properties['wf_increment'])
            t_start = float(data_object.properties['wf_start_offset'])
            t_end   = float(data_object.properties['wf_samples']) * t_step + t_start - t_step/2
            t_list  = np.arange(t_start, t_end, t_step)

            data_dict[par] = [t_list, np.array(data_object)]

    return data_dict


def calibration():
    param_list = ['LC', 'PS']
    cal_original = read_tdms('test_data_2021/Calibration_211217_105057.tdms', param_list)

    res = 1000
    delta = 20
    max_diffs = [0.3 * 10**(-5), 0.3 * 10**(-4)]
    ranges_list = [[[-0.00007, 0.00001], [-0.0001, 0.00001], [-0.00021, 0.00002]], [[0.004, 0.0001], [0.0049, 0.0001]]]

    values = [[0, 3.938 * 9.80665, 17.593 * 9.80665], [0, 7 * 10**5]]  # N, Pa
    slope_list = []
    offset_list = []

    for param, max_diff, ranges, val in zip(param_list, max_diffs, ranges_list, values):
        cal_t    = cal_original[param][0][::res]
        cal_data = cal_original[param][1][::res]

        diff = abs(cal_data[:-delta] - cal_data[delta:])
        filtered = cal_data[np.where(diff < max_diff)]

        avg_list = [np.average(filtered[np.where(abs(filtered - r[0]) < r[1])]) for r in ranges]

        # plt.plot(cal_t, cal_data)
        # # plt.scatter(cal_t[np.where(diff < max_diff)], filtered, s=3)
        # for avg in avg_list:
        #     plt.hlines(avg, cal_t[0], cal_t[-1], colors='r')
        # plt.show()

        slope1 = (val[1]  - val[0])  / (avg_list[1]  - avg_list[0])
        slope2 = (val[-1] - val[0])  / (avg_list[-1] - avg_list[0])
        slope3 = (val[-1] - val[-2]) / (avg_list[-1] - avg_list[-2])
        slope_list.append((slope1 + slope2 + slope3) / 3)

        offset_list.append(avg_list[0])

    return slope_list, offset_list


def convert(filename, begin):
    params = ['LC', 'PS']
    filtered_dict = {}

    slopes, offsets = calibration()
    data = read_tdms(filename, params)

    for param, s, off in zip(params, slopes, offsets):
        t_range = np.where(np.logical_and(data[param][0] > begin, data[param][0] < begin + 15))
        data_t = data[param][0][t_range][::100]
        data_y = data[param][1][t_range][::100] - off

        filtered = data[param][1][t_range][::100] - off
        n_points = 30
        if param == 'LC':
            n_points = 1

        for i in range(n_points):
            filtered[:-n_points] += data_y[i:-n_points + i]
        filtered /= (n_points + 1)

        if param == 'LC':
            dt = data_t[1] - data_t[0]
            impulse = np.zeros(len(filtered[:-n_points]) + 1)

            for j in range(1, len(impulse)):
                impulse[j] = impulse[j-1] + filtered[j-1] * dt
            impulse = impulse[1:]

            filtered_dict['IM'] = [data_t[:-n_points], impulse * s]
            plt.plot(*filtered_dict['IM'])
            plt.show()

        filtered_dict[param] = [data_t[:-n_points], filtered[:-n_points] * s]
        plt.plot(*filtered_dict[param])
        plt.show()

    return filtered_dict


def save_data(filename_list, t_begin_list):
    with open('Filtered_Data.txt', 'w') as file:
        for filename, t_begin in zip(filename_list, t_begin_list):
            file.write(filename + '\n')
            filtered_data = convert(filename, t_begin)

            for p in filtered_data.keys():
                file.write(p + '_time: ' + str(filtered_data[p][0].tolist())[1:-1] + '\n')
                file.write(p + ': '      + str(filtered_data[p][1].tolist())[1:-1] + '\n')


def read_data(filename):
    with open('Filtered_Data.txt') as file:
        start = np.infty
        data_dict = {}

        for i, line in enumerate(file):
            if filename in line.strip():
                start = i + 1

            if i >= start + 6:
                break

            if i >= start:
                param, data = line.strip().split(': ')
                data_dict[param] = np.array(data.split(', '), dtype=float)

    return data_dict


if __name__ == '__main__':
    # convert('test_data_2021/Config2_211221_132537.tdms', 34)
    # convert('test_data_2021/ReferenceMotor_211221_114135.tdms', 222)

    # save_data(['test_data_2021/Config2_211221_132537.tdms', 'test_data_2021/ReferenceMotor_211221_114135.tdms'],
    #           [34, 222])

    data_ = read_data('Con')
    print(data_.keys())
    print(data_['PS'][3])
