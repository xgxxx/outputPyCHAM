''' Module to plot SOA temporal trends correlated with solar insolation data '''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def SOA_insolation(time, time_interval, file_name, sheet, start_time, particulate_phase_mass):
    ##### get needed length
    length = int(time[-1] / time_interval)

    ### get solar information (remove the data at t=0)
    insolation = pd.ExcelFile(file_name)
    insolation = insolation.parse(sheet)
    insolation = insolation.to_numpy()
    insolation_need = insolation[insolation[:, 3] == pd.Timestamp(start_time)]
    ref = np.array([True, True, True, True, True, True, True, True, True, True, True, True, True, False])
    ref2 = (insolation == insolation_need[0])
    start = np.where((ref2 == ref).all(axis=1))[0][0]
    solar = np.array([insolation[start + i + 1][11] for i in range(length)])

    ###calculate total SOA in all size bins
    SOA_total = np.sum(particulate_phase_mass, axis=0)
    SOA_total_interval = np.zeros(length)
    for i in range(length):
        SOA_total_interval[i] = SOA_total[time_interval * (i+1)]

    ###scatter plot and save
    plt.figure(1)
    plt.scatter(solar, SOA_total_interval)
    plt.xlabel("Solar Radiation (W/m2)")
    plt.ylabel("SOA concentration (ug/m3)")
    plt.title("SOA temporal trends (time interval: " + str(time_interval) + " min)")
    plt.savefig("SOA temporal trends-Solar radiation")

