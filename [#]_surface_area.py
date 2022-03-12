''' Module to calculate [#](number concentration) and surface area'''

import numpy as np
import pandas as pd
from math import pi


def number_concentration():
    return particulate_phase

def surface_area(radius, components_number, particulate_phase):
    ### surface area
    particulate_phase_need = particulate_phase[:, 2:particulate_phase.shape[1]]
    particulate_phase_need = particulate_phase_need.astype(np.float)
    surface = np.zeros((len(particulate_phase_need), len(particulate_phase_need[0])))
    for i in range(len(particulate_phase_need)):
        for j in range(len(particulate_phase_need[0])):
            surface[i][j] = (particulate_phase_need[i][j]) * pi * (radius[int(i / components_number)] * 2) ** 2
    pd.DataFrame(surface).to_csv('surface area.csv')


