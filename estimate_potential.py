''' estimate potential of formation for each gas '''

import numpy as np
import pandas as pd

def estimate_potential(gas_phase,total,components_number,species):
    ### gas phase
    particulate_phase_gas = gas_phase
    particulate_phase_gas = particulate_phase_gas[:, 2:particulate_phase_gas.shape[1]]
    particulate_phase_gas = particulate_phase_gas.astype(np.float)

    ### particle phase
    particulate_phase_particular = total

    ### components_number & time
    times = len(total[0])

    ###potential
    potential = np.zeros((components_number, times))
    for i in range(components_number):
        for j in range(times):
            if particulate_phase_gas[i][j] == 0:
                potential[i][j] = 'nan' ###gas phase = 0
            else:
                potential[i][j] = particulate_phase_particular[i][j] / particulate_phase_gas[i][j]
    potential = pd.concat([pd.DataFrame(species),
                           pd.DataFrame(potential)],
                          axis=1)
    potential.to_csv('GHG.csv')
