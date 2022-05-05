''' Module to plot compare simulation result and observation result  '''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def compare(skip_num, SOA_list, start_time, end_time, time_column, particulate_phase_mass_SOA, time_interval, file_name, sheet):
    # start_time: start of needed time, this time will be included in calculation \
    # the format should match with that in time_column, given by user (e.g. in time column "Date&Time", \
    # the format of time is "year-month-date hour:minute:second", so one example of start time is "2020-11-19 17:05:00")
    # end_time: end of needed time, this time will be included in calculation \
    # the format should match with that in time_column, given by user (e.g. in time column "Date&Time", \
    # the format of time is "year-month-date hour:minute:second", so one example of end time is "2020-11-19 18:05:00")
    # particulate_phase_mass_SOA: contains partitulate phase mass concentration (unit: ug/m3)\
    # for each SOA component in each size bin and each time, given in post processing file
    # time_interval: time interval of observation data (unit: minute),  given by user
    # file_name: name of the file that contains observation data, given by user (e.g. "5-min SOA.xlsx")
    # sheet: name of the sheet that contains observation data in the above file, given by user (e.g. "AMS")

    ### extract obervation data
    observation = pd.ExcelFile(file_name)
    observation = observation.parse(sheet, skiprows=skip_num)
    observation = observation[:3100]
    observation = observation.loc[observation[time_column] >= start_time]
    observation = observation.loc[observation[time_column] <= end_time]
    observation_time = observation[time_column]
    observation_data = []
    for i in SOA_list:
        observation_data = observation_data+observation[i]

    ### simulation data
    SOA_total = particulate_phase_mass_SOA[1:, 2:]
    SOA_total = SOA_total.astype(np.float)
    SOA_total = np.sum(SOA_total, axis=0)
    SOA_total_average = [SOA_total[0]]

    for i in range(len(observation)-1):
        SOA_total_average.append(np.mean(SOA_total[1+i*time_interval:1+(i+1)*time_interval]))

    ###convert to csv and plot
    comparison = pd.DataFrame({'observation time': observation_time,
                               'simulation (ug/m3)': SOA_total_average,
                               'observation (ug/m3)': observation_data})
    comparison.to_csv("comparison between simulation and observation1.csv")
