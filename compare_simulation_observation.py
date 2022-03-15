''' Module to plot compare simulation result and observation result  '''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def compare(start_end, particulate_phase_mass, time, time_interval, file_name, observation_start_time):
    if start_end==0:
        ##### get needed length
        length = int(time[-1] / time_interval)
    else:
        ##### get needed length
        length = int((time[-1]-5) / time_interval)
    ###calculate total SOA in all size bins (average, time_interval)
    SOA_total = np.sum(particulate_phase_mass, axis=0)
    SOA_total_average = np.zeros(length)
    for i in range(length):
        SOA_total_average[i] = np.mean(SOA_total[(start_end+1) + 10 * i:(start_end+1) + 10 * (i + 1)])

    ### extract from observation data
    observation = pd.read_excel(file_name)
    observation = observation.to_numpy()
    observation = observation[2:, :9]
    observation_need = observation[observation[:, 0] == observation_start_time]
    observation_need = observation_need[0]
    start = np.where((observation == observation_need).all(axis=1))[0][0]
    SOA_observation_need = np.array([observation[start + i][3] + observation[start + i][5] + observation[start + i][7] for i in range(length)])
    observation_time = np.array([observation[start + i][1] for i in range(length)])

    ###convert to csv and plot
    comparison = pd.DataFrame({'simulation': SOA_total_average,
                               'observation': SOA_observation_need})
    comparison.to_csv("comparison between simulation and observation.csv")

    plt.figure(1)
    plt.plot(observation_time, SOA_total_average, 'r')  # red line is for simulation result
    plt.plot(observation_time, SOA_observation_need, 'g')  # green line is for observation result
    plt.title("comparison of simulation and observation(red:simulation; green:observation)")
    plt.xlabel("time (min)")
    plt.ylabel("SOA concentration (ug/m3)")
    plt.savefig("comparison (simulation and observation)")

