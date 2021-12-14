import numpy as np


def read_datafile():
    params_dict = {}
    data_dict = {}
    data = False
    i_next_data = 1000

    with open('TRP_exercise.txt') as file:
        for i, line in enumerate(file):
            if i <= 2:
                key, value = line.strip().split(': ')
                params = value.split(' "],[" ')

                params[0] = params[0].strip('[" ')
                params[-1] = params[-1].strip(' "]')

                params_dict[key] = params

            if i > 2:
                if '###' in line:
                    if data:
                        data_dict[data['Pressure']] = data

                    i_next_data = i + 1
                    data = {}

                if 0 <= i - i_next_data < 26:
                    data[params_dict['States'][i - i_next_data]] = np.float(line.strip())

                if 26 + 1 <= i - i_next_data < 26 + 1 + 5:
                    data[params_dict['Optimum'][i - i_next_data - (26 + 1)] + '_opt'] = np.float(line.strip())

                if 26 + 1 + 5 + 1 <= i - i_next_data < 26 + 1 + 5 + 1 + 4:
                    data[params_dict['Vacuum'][i - i_next_data - (26 + 1 + 5 + 1)] + '_vac'] = np.float(line.strip())

    return params_dict, data_dict


def reform(data_dict):
    new_data_dict = data_dict[1.5]

    for param in new_data_dict.keys():
        new_data_dict[param] = [data[param] for data in data_dict.values()]

    return new_data_dict


def linearize(param_str, p):
    data_dict = reform(read_datafile()[1])

    p /= 1000000
    p_lst = data_dict['Pressure']

    i = 0
    while p_lst[i] < p:
        i += 1

    dif = p - p_lst[i-1]
    step = p_lst[i] - p_lst[i-1]

    smaller = data_dict[param_str][i-1]
    bigger  = data_dict[param_str][i]

    return smaller + (bigger - smaller) * (dif / step)


if __name__ == '__main__':
    # print(read_datafile()[0])
    # print(read_datafile()[1])

    new_data = reform(read_datafile()[1])
    print(new_data)

    print(linearize('Temperature', 1.0 * 10**6))
