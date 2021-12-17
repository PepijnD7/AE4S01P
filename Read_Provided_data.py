import numpy as np
import time
import matplotlib.pyplot as plt


def read_given_datafile():
    t_start = time.time()
    given_data = {}
    params_lst = []

    with open('given_test_data.txt') as file:
        for i, line in enumerate(file):
            if line.strip() == '':
                start_line = i

            if i - start_line == 1:
                param = line.strip()
                params_lst.append(param)
                given_data[param] = []

            if i - start_line > 1:
                given_data[param].append(float(line.strip()))

            t = time.time()
            if t - t_start > 10:
                print(param, len(given_data[param]))
                break

    for key in params_lst:
        given_data[key] = np.array(given_data[key])
        # print(len(given_data[key]))

    new_param = 'Time [s]'
    given_data[new_param] = given_data[params_lst[1]] / 1024000
    given_data[new_param] -= given_data[new_param][0] + 37.24
    params_lst[1] = new_param

    new_param = 'Pressure [Pa]'
    given_data[new_param] = (given_data[params_lst[0]] * 6250 - 25) * 10**5
    given_data[new_param] -= given_data[new_param][0]
    params_lst[0] = new_param

    return given_data, params_lst


def plot_given_data(our_data=None):
    data, params = read_given_datafile()
    r = [36000, 44000]  # data range

    plt.plot(data[params[1]][r[0]:r[1]], data[params[0]][r[0]:r[1]])
    plt.xlabel(params[1])
    plt.ylabel(params[0])
    plt.grid()
    plt.show()

    if our_data:
        t, p = our_data

        plt.plot(t, p, label='Simulation')
        plt.plot(data[params[1]][r[0]:r[1]], data[params[0]][r[0]:r[1]], label='Test Data')
        plt.xlabel(params[1])
        plt.ylabel(params[0])
        plt.legend()
        plt.grid()
        plt.show()


if __name__ == '__main__':
    plot_given_data()
