''' Module to calculate [#](number concentration) and surface area'''

import numpy as np
import pandas as pd
from math import pi


def surface_for_each_component(radius, components_number, particulate_phase):
    # radius: radius of each size bin at each time, extracted from PyCHAM output size_bin_radius file, \
    # given in post processing file
    # components_number: number of components, given in post processing file
    # particulate_phase: includes mass concentration corresponding to each component's name \
    # in each size bin at each time, given in post processing file

    surface = particulate_phase
    for component in range(len(particulate_phase)):
        for times in range(len(particulate_phase[0])):
            # if times = 0 or 1, no change is needed as they indicate the component's name and size bin number
            # only do changes when times >= 2
            if times >= 2:
                # calculate by given formula, need to transpose radius, and multiply radius by 2 to convert to diameter
                surface[component][times] = (float(particulate_phase[component][times])) * pi * (
                            radius.transpose()[int(component / components_number)][times-2] * 2) ** 2

    # col_name: name of each column
    col_name = np.array([[' ', ' ']+[str(i)+' minute' for i in range(len(time))]])
    surface = np.vstack((col_name, surface))
    pd.DataFrame(surface).to_csv('surface_for_each_component.csv')


def surface_for_total(radius, bin_number, components_number_SOA, particulate_phase_SOA):
    # radius: radius of each size bin at each time, extracted from PyCHAM output size_bin_radius file, \
    # given in post processing file
    # components_number_SOA: number of SOA components, given in post processing file
    # particulate_phase_mass_SOA: contains partitulate phase mass concentration\
    # for each SOA component in each size bin and each time, given in post processing file

    # exclude first 2 columns which contains component's name and size bin number
    pp_SOA = particulate_phase_SOA[:, 2:particulate_phase_SOA.shape[1]]
    pp_SOA = pp_SOA.astype(np.float)

    col_name = [' '] + [str(i) + ' minute' for i in range(len(time))]
    surface = [col_name]

    for bins in range(bin_number):
        # calculate the total SOA mass concentration in each size bin at each time
        temp_sum = pp_SOA[components_number_SOA * bins:components_number_SOA * (bins + 1)]
        temp_sum = np.sum(temp_sum, axis=0)

        # add 1 column which indicates the size bin number
        bin_surface = ['p' + str(bins+1)]

        for times in range(len(time)):
            bin_surface.append((temp_sum[times]) * pi * (
                            radius.transpose()[bins][times] * 2) ** 2)
        surface.append(bin_surface)
    pd.DataFrame(surface).to_csv('surface_for_total.csv')
