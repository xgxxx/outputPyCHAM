''' Module to plot compare simulation result and observation result  '''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def compare(smaller, length, start_time, end_time, time_column, particulate_phase_mass_SOA, time_interval, file_name, sheet):
    # smaller: check if simulation start time is smaller than start time of observation, \
    # if it is, smaller = 1, else smaller = 0
    # length: number of data that need to be compared (due to the input file format)
    # start_time: start time of observation, given by user
    # end_time: end time of observation, given by user
    # particulate_phase_mass_SOA: contains partitulate phase mass concentration\
    # for each SOA component in each size bin and each time, given in post processing file
    # time_interval: time interval of observation data,  given by user
    # file_name: name of the file that contains observation data, given by user
    # sheet: name of the sheet that contains observation data in the above file, given by user

    ### extract obervation data
    observation = pd.ExcelFile(file_name)
    observation = observation.parse(sheet, skiprows=2)
    observation = observation.loc[observation[time_column] >= start_time]
    observation = observation.loc[observation[time_column] <= end_time]
    observation = observation[:length]
    observation_time = observation[time_column]
    observation_data = observation["MOOOA"] + observation["LOOOA"] + observation["OOA"]
    SOA_total = np.sum(particulate_phase_mass_SOA, axis=0)
    SOA_total_average = []

    if smaller == 1:
        start = 5 - int(start_time[-5:-3])+1
        SOA_total_average.append(np.mean(SOA_total[:start]))
    else:
        start = 1
        SOA_total_average.append(' ')
        
    for i in range(length-1):
        SOA_total_average.append(np.mean(SOA_total[start+i*time_interval:start+(i+1)*time_interval]))

    ###convert to csv and plot
    comparison = pd.DataFrame({'observation time': observation_time,
                               'simulation (ug/m3)': SOA_total_average,
                               'observation (ug/m3)': observation_data})
    comparison.to_csv("comparison between simulation and observation1.csv")
