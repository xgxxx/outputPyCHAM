''' estimate potential of formation for each gas '''

import numpy as np
import pandas as pd

def estimate_potential(gas_phase_mass_SOA,total_particulate_mass_species_time_SOA):
    # gas_phase_mass_SOA: mass concentration for each component at each time (ug/m3), given in post processing file \
    # (e.g. 6.1203140008331e-07 ug/m3)
    # total_particulate_mass_species_time_SOA: total mass concentration for each component in all size bins \
    # at each time (ug/m3), given in post processing file (5.22207160752873e-18 ug/m3)

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
    potential = total_particulate_mass_species_time_SOA[1:]
    for i in range(len(potential)):
        potential[i] = np.concatenate((potential[i][:1], potential_SOA[i]))
    col_name = np.array([[' '] + [str(i) + ' minute' for i in range(len(total_particulate_mass_species_time_SOA[0])-1)]])
    potential = np.vstack((col_name, potential))
    pd.DataFrame(potential).to_csv('SOA_potential.csv')
