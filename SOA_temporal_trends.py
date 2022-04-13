''' Module to plot SOA temporal trends correlated with solar insolation data '''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def SOA_insolation(time_interval, file_name, sheet, time_column, start_time, end_time, insolation_column, particulate_phase_mass_SOA):
    # time_interval: time interval of observed data, given by user
    # file_name: name of the file that contains insolation data, given by user
    # sheet: name of the sheet that contains insolation data in the above file, given by user
    # time_column: name of the column that includes detailed observation time, given by user
    # start_time: start of needed time, given by user
    # end_time: end of needed time, given by user
    # insolation_column: name of the column that includes insolation data, given by user
    # particulate_phase_mass_SOA: given in main file, contains partitulate phase mass concentration\
    # for each SOA component in each size bin and each time


    ### get insolation information (remove the data at t=0)
    insolation = pd.ExcelFile(file_name)
    insolation = insolation.parse(sheet)
    insolation = insolation.loc[insolation[time_column] > start_time]
    insolation = insolation.loc[insolation[time_column] <= end_time]
    insolation = insolation[insolation_column]
    insolation = insolation.to_numpy()


    ###calculate total SOA in all size bins
    SOA_total = np.sum(particulate_phase_mass_SOA, axis=0)
    SOA_total_interval = np.zeros(len(insolation))
    for i in range(len(insolation)):
        SOA_total_interval[i] = np.mean(SOA_total[1+time_interval*i:1+time_interval*(i+1)])
    ###scatter plot and save
    plt.figure(1)
    plt.scatter(insolation, SOA_total_interval)
    plt.xlabel("Solar Radiation (W/m\u00b2)")
    plt.ylabel("SOA concentration (\u03BCg/m\u00b3)")
    plt.title("SOA temporal trends (time interval: " + str(time_interval) + " minutes)")
    plt.show()
    plt.savefig("SOA temporal trends-Solar radiation")

