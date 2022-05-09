''' Module to plot SOA temporal trends correlated with solar insolation data '''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def SOA_insolation(time_interval, file_name, sheet, time_column, start_time, end_time, insolation_column, particulate_phase_mass_SOA):
    # time_interval: time interval of observed data (unit: minute), given by user. It should be match to simulation update step, \
    # for example, if time interval is 5 minute and simulation update step is 7 minute, they are NOT match
    # file_name: name of the file that contains insolation data, given by user (e.g. "Meteorology_2020.xlsx")
    # sheet: name of the sheet that contains insolation data in the above file, given by user \
    # (e.g. "Our weather station")
    # time_column: name of the column that contains detailed observation time, given by user (e.g. "Date&Time")
    # start_time: start of needed time, this time will be included in calculation \
    # the format should match with that in time_column, given by user (e.g. in time column "Date&Time", \
    # the format of time is "year-month-date hour:minute:second", so one example of start time is "2020-11-19 17:00:00")
    # end_time: end of needed time, this time will be included in calculation \
    # the format should match with that in time_column, given by user (e.g. in time column "Date&Time", \
    # the format of time is "year-month-date hour:minute:second", so one example of end time is "2020-11-19 17:55:00")
    # insolation_column: name of the column that contains insolation data, given by user (e.g. "Solar Radiation (W/m2)")
    # particulate_phase_mass_SOA: contains particulate phase mass concentration for each SOA component \
    # in each size bin and each time (ug/cm3), given in post processing file (e.g. 6.752614986953e-26 ug/cm3)


    ### get insolation information
    insolation = pd.ExcelFile(file_name)
    insolation = insolation.parse(sheet)
    insolation = insolation.loc[insolation[time_column] >= start_time]   # start time is included in calculation)
    insolation = insolation.loc[insolation[time_column] <= end_time]   # end time is included in calculation)
    insolation_time = insolation[time_column]
    insolation_time = list(insolation_time)
    insolation = insolation[insolation_column]
    insolation = insolation.to_numpy()


    ### calculate total SOA in all size bins
    SOA_total = particulate_phase_mass_SOA[1:, 2:]
    SOA_total = SOA_total.astype(np.float)
    SOA_total = np.sum(SOA_total, axis=0)
    num_data = int(len(SOA_total)/time_interval)
    SOA_total_interval = np.zeros(num_data)
    for i in range(num_data):
        SOA_total_interval[i] = np.mean(SOA_total[time_interval*i:time_interval*(i+1)])

    ### scatter plot and save
    plt.figure(1)
    plt.scatter(insolation, SOA_total_interval)
    plt.xlabel("Solar Radiation (W/m\u00b2)")
    plt.ylabel("SOA concentration (\u03BCg/m\u00b3)")
    plt.title("SOA temporal trends (time interval: " + str(time_interval) + " minutes)")
    plt.savefig("SOA temporal trends-Solar radiation")

    ### convert to csv file
    SOA_temporal_trends = [['Time', 'SOA (ug/m3)', 'Insolation (W/m2)']]
    for i in range(len(insolation)):
        SOA_temporal_trends.append([insolation_time[i], SOA_total_interval[i], insolation[i]])
    pd.DataFrame(SOA_temporal_trends).to_csv("SOA_temporal_trends.csv", index=False, header=False)
