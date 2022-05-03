''' estimate potential of formation for each gas '''

import numpy as np
import pandas as pd

def estimate_potential(gas_phase_mass_SOA,total_particulate_mass_species_time_SOA):
    ### gas phase
    gas = gas_phase_mass_SOA
    gas = gas[1:, 2:]
    gas = gas.astype(np.float)

    ### particle phase
    particulate = total_particulate_mass_species_time_SOA
    particulate = particulate[1:, 1:]
    particulate = particulate.astype(np.float)

    ###potential
    potential_SOA = particulate / gas
    potential = gas_phase_mass_SOA[1:, 1:]
    for i in range(len(potential)):
        potential[i] = np.concatenate((potential[i][:1], potential_SOA[i]))
    col_name = np.array([[' '] + [str(i) + ' minute' for i in range(len(time))]])
    potential = np.vstack((col_name, potential))
    pd.DataFrame(potential).to_csv('potential.csv')
