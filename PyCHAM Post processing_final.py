#####################################
## BN ##
import numpy as np
import pandas as pd 

# Set working directory - input folder of interest
import os
path = r'C:\Users\janet\Documents\nBox\0. SOAFP-PyCHAM output\PR3_simulation\Simulation results\20Nov_1700_1759_89VOCs_20bins'
os.chdir(path)

##### Size bin radius (unit: um & mm)
radius_file_path = r"size_bin_radius"
radius = open(radius_file_path, "r+")
radius = radius.readlines()
radius_um = pd.DataFrame([radius[i].split(",") for i in range(2, len(radius))])
radius_um.to_csv('Radius (\u03BCm).csv')
radius = [radius[i].split(",") for i in range(2, len(radius))]
radius = np.array(radius).astype(float)
radius_mm = radius*10**(-3)
pd.DataFrame(radius_mm).to_csv('Radius (mm).csv')

##### Size bin diameter (unit: um & mm)
diameter = radius*2
diameter_um = diameter
pd.DataFrame(diameter_um).to_csv('Diameter (\u03BCm).csv')
diameter_mm = diameter*10**(-3)
pd.DataFrame(diameter_mm).to_csv('Diameter (mm).csv')

##### time (unit: min)
time_path = r"time"
time = open(time_path, "r+")
time = time.readlines()
time = [time[i] for i in range(1, len(time))]
time = np.array([float(i)/60 for i in time])

##### Concentration information
# Input filepath here #
in_file_path = r"concentrations_all_components_all_times_gas_particle_wall"
information_file_path = r"model_and_component_constants"
chamber_environment_file_path = r"chamber_environmental_conditions"

# Read concentration data
concentration_file = open(in_file_path, "r+")
concentration_file = concentration_file.readlines()
species = np.array([x.split("_") for x in concentration_file[1].split(",")])
species_1 = [species[x][0] for x in range(len(species))]
species_2 = [species[x][1] for x in range(len(species))]
species = np.vstack((species_1,species_2)).transpose()
concentration = np.array([concentration_file[i].split(",") for i in range(2, len(concentration_file))]).transpose()
concentration = np.concatenate((species,concentration), axis=1)

# Extract information for each species
# Molecular weight (g/mol)
information = open(information_file_path, "r+")
information = information.readlines()
molecular_weight = information[2].split(",")
# Fix the first and last
molecular_weight = molecular_weight[1:len(molecular_weight)] 
molecular_weight[0] = molecular_weight[0][1:len(molecular_weight[0])]
molecular_weight[-1] = molecular_weight[-1][0:len(molecular_weight[-1])-2]
# Convert to float
molecular_weight = [float(i) for i in molecular_weight]

factor = information[8].split(",")
factor = factor[1:len(factor)]
factor = [x.strip(" []\n") for x in factor]
factor = np.array(factor).astype(float).transpose()

## Extract radicals
# Peroxy Radical
organic_peroxy_radical_index = information[4].split(",")
organic_peroxy_radical_index = organic_peroxy_radical_index[1:len(organic_peroxy_radical_index)]
organic_peroxy_radical_index[0] = organic_peroxy_radical_index[0][1:len(organic_peroxy_radical_index[0])]
organic_peroxy_radical_index[-1] = organic_peroxy_radical_index[-1][0:len(organic_peroxy_radical_index[-1])-2]
# Convert to int
organic_peroxy_radical_index = [int(i) for i in organic_peroxy_radical_index]

# Alkoxy Radical
organic_alkoxy_radical_index = information[5].split(",")
organic_alkoxy_radical_index = organic_alkoxy_radical_index[1:len(organic_alkoxy_radical_index)]
organic_alkoxy_radical_index[0] = organic_alkoxy_radical_index[0][1:len(organic_alkoxy_radical_index[0])]
organic_alkoxy_radical_index[-1] = organic_alkoxy_radical_index[-1][0:len(organic_alkoxy_radical_index[-1])-2]
# Convert to int
organic_alkoxy_radical_index = [int(i) for i in organic_alkoxy_radical_index]

#Extract size_bin_number and components_number
bin_number = int(information[0].split(",")[1])
components_number = int(information[1].split(",")[1])

##### non SOA components
nonSOA = ['core', 'H2O', 'AMMSUL']
components_number_SOA = components_number-len(nonSOA)

#Extract saturation vapor pressure at 298.15K of species
saturation_vapor_pressure = information[13].split(",")
saturation_vapor_pressure = saturation_vapor_pressure[1:len(saturation_vapor_pressure)]
saturation_vapor_pressure[0] = saturation_vapor_pressure[0][1:len(saturation_vapor_pressure[0])]
saturation_vapor_pressure[-1] = saturation_vapor_pressure[-1][0:len(saturation_vapor_pressure[-1])-2]
# Convert to int
saturation_vapor_pressure = [float(i) for i in saturation_vapor_pressure]

# Extract gas phase concentration
gas_phase = concentration[concentration[:,1] == "g"]
gas_phase[0,0] = gas_phase[0,0][2:]

# species corresponding to molecular weight
species = gas_phase[0:(len(gas_phase)),0].tolist()

# Read chamber environment temperature
environment = open(chamber_environment_file_path, "r+")
environment = environment.readlines()
environment = environment[1:len(environment)]
temperature = [environment[i].split(",")[0] for i in range(len(environment))]
temperature = np.array(temperature).transpose().astype(float)
rh = [environment[i].split(",")[2] for i in range(len(environment))]
rh = np.array(rh).transpose().astype(float)

# Extract particulate phase concentration (molecules/cm3)
particulate_phase = concentration[concentration[:,1] != "g"]
particulate_phase[-1][1] = particulate_phase[-1][1][:-1]

particulate_phase_all = particulate_phase
col_name = np.array([[' ', ' ']+[str(i)+' minute' for i in range(len(time))]])
particulate_phase_all = np.vstack((col_name, particulate_phase_all))

particulate_phase_SOA = particulate_phase_all
for i in nonSOA:
    particulate_phase_SOA = particulate_phase_SOA[particulate_phase_SOA[:, 0] != ' '+i]

particulate_phase_nonSOA = np.array([i for i in particulate_phase_all if i[0][1:] in nonSOA])
col_name = np.array([[' ', ' ']+[str(i)+' minute' for i in range(len(time))]])
particulate_phase_nonSOA = np.vstack((col_name, particulate_phase_nonSOA))

pd.DataFrame(particulate_phase_all).to_csv("particulate_phase_all (molecules.(cc (air))\u207B\u00b9).csv")
pd.DataFrame(particulate_phase_SOA).to_csv("particulate_phase_SOA (molecules.(cc (air))\u207B\u00b9).csv")
pd.DataFrame(particulate_phase_nonSOA).to_csv("particulate_phase_nonSOA (molecules.(cc (air))\u207B\u00b9).csv")

# Convert from molecules/cm3 to ppb
particulate_phase_ppb = particulate_phase[:, 2:particulate_phase.shape[1]]
particulate_phase_ppb = particulate_phase_ppb.astype(np.float)
particulate_phase_ppb = particulate_phase_ppb / factor

particulate_phase_ppb_all = particulate_phase
for i in range(len(particulate_phase)):
    particulate_phase_ppb_all[i] = np.concatenate((particulate_phase[i][:2], particulate_phase_ppb[i]))
col_name = np.array([[' ', ' ']+[str(i)+' minute' for i in range(len(time))]])
particulate_phase_ppb_all = np.vstack((col_name, particulate_phase_ppb_all))

particulate_phase_ppb_SOA = particulate_phase_ppb_all
for i in nonSOA:
    particulate_phase_ppb_SOA = particulate_phase_ppb_SOA[particulate_phase_ppb_SOA[:, 0] != ' '+i]

particulate_phase_ppb_nonSOA = np.array([i for i in particulate_phase_ppb_all if i[0][1:] in nonSOA])
col_name = np.array([[' ', ' ']+[str(i)+' minute' for i in range(len(time))]])
particulate_phase_ppb_nonSOA = np.vstack((col_name, particulate_phase_ppb_nonSOA))

pd.DataFrame(particulate_phase_ppb_all).to_csv("particulate_phase_all (ppb).csv")
pd.DataFrame(particulate_phase_ppb_SOA).to_csv("particulate_phase_SOA (ppb).csv")
pd.DataFrame(particulate_phase_ppb_nonSOA).to_csv("particulate_phase_nonSOA (ppb).csv")

# repeat molecular weight array by no. of times of size bins
temp = np.tile(molecular_weight, int(len(particulate_phase_ppb) / len(molecular_weight)))
# Convert from ppb to ug/m3
particulate_phase_mass = particulate_phase_ppb * 12.187 / temperature
particulate_phase_mass = (particulate_phase_mass.transpose() * temp.transpose()).transpose()

particulate_phase_mass_all = particulate_phase
for i in range(len(particulate_phase)):
    particulate_phase_mass_all[i] = np.concatenate((particulate_phase[i][:2], particulate_phase_mass[i]))
col_name = np.array([[' ', ' ']+[str(i)+' minute' for i in range(len(time))]])
particulate_phase_mass_all = np.vstack((col_name, particulate_phase_mass_all))

particulate_phase_mass_SOA = particulate_phase_mass_all
for i in nonSOA:
    particulate_phase_mass_SOA = particulate_phase_mass_SOA[particulate_phase_mass_SOA[:, 0] != ' '+i]

particulate_phase_mass_nonSOA = np.array([i for i in particulate_phase_mass_all if i[0][1:] in nonSOA])
col_name = np.array([[' ', ' ']+[str(i)+' minute' for i in range(len(time))]])
particulate_phase_mass_nonSOA = np.vstack((col_name, particulate_phase_mass_nonSOA))

pd.DataFrame(particulate_phase_mass_all).to_csv("particulate_phase_all (\u03BCg.m\u207B\u00b3).csv")
pd.DataFrame(particulate_phase_mass_SOA).to_csv("particulate_phase_SOA (\u03BCg.m\u207B\u00b3).csv")
pd.DataFrame(particulate_phase_mass_nonSOA).to_csv("particulate_phase_nonSOA (\u03BCg.m\u207B\u00b3).csv")

### total mass concentration (ug/m3) per species per size bin
total_particulate_mass_species_sizebin = np.sum(particulate_phase_mass, axis=1)
total_particulate_mass_species_sizebin = total_particulate_mass_species_sizebin.reshape(bin_number,components_number).transpose()
# total mass concentration (ug/m3) per species per size bin for all components
total_particulate_mass_species_sizebin_all = particulate_phase[:components_number]
total_particulate_mass_species_sizebin_all = total_particulate_mass_species_sizebin_all[:,:bin_number+1]
for i in range(components_number):
    total_particulate_mass_species_sizebin_all[i] = np.concatenate((particulate_phase[i][:1], total_particulate_mass_species_sizebin[i]))
col_name = np.array([['Species']+['p '+str(i+1) for i in range(bin_number)]])
total_particulate_mass_species_sizebin_all = np.vstack((col_name, total_particulate_mass_species_sizebin_all))
pd.DataFrame(total_particulate_mass_species_sizebin_all).to_csv('perSpecies_perSizebin_all (\u03BCg.m\u207B\u00b3).csv')

# total mass concentration (ug/m3) per species per size bin for SOA components
total_particulate_mass_species_sizebin_SOA = total_particulate_mass_species_sizebin_all
for i in nonSOA:
    total_particulate_mass_species_sizebin_SOA = total_particulate_mass_species_sizebin_SOA[total_particulate_mass_species_sizebin_SOA[:, 0] != ' '+i]
pd.DataFrame(total_particulate_mass_species_sizebin_SOA).to_csv('perSpecies_perSizebin_SOA (\u03BCg.m\u207B\u00b3).csv')

# total mass concentration (ug/m3) per species per size bin for nonSOA components
total_particulate_mass_species_sizebin_nonSOA = np.array([i for i in total_particulate_mass_species_sizebin_all if i[0][1:] in nonSOA])
col_name = np.array([['Species']+['p '+str(i+1) for i in range(bin_number)]])
total_particulate_mass_species_sizebin_nonSOA = np.vstack((col_name, total_particulate_mass_species_sizebin_nonSOA))
pd.DataFrame(total_particulate_mass_species_sizebin_nonSOA).to_csv('perSpecies_perSizebin_nonSOA (\u03BCg.m\u207B\u00b3).csv')

### total mass concentration (ug/m3) per species
total_particulate_mass_species = np.sum(total_particulate_mass_species_sizebin, axis = 1)
# total mass concentration (ug/m3) per species per time
total_particulate_mass_species_time = np.array_split(particulate_phase_mass, int(len(particulate_phase_ppb)/len(molecular_weight)))
total = np.zeros(shape=(len(molecular_weight),len(environment)))
for items in total_particulate_mass_species_time:
    total = total + items

##### particle_number_concentration_dry_number
particle_number_concentration_dry_file_path = r"particle_number_concentration_dry"
particle_number_concentration_dry = open(particle_number_concentration_dry_file_path, "r+")
particle_number_concentration_dry = particle_number_concentration_dry.readlines()
particle_number_concentration_dry = pd.DataFrame([particle_number_concentration_dry[i].split(",") for i in range(2, len(particle_number_concentration_dry))])
particle_number_concentration_dry.to_csv('particle_number_concentration_dry_number.csv')

##### particle_number_concentration_wet_number
particle_number_concentration_wet_file_path = r"particle_number_concentration_wet"
particle_number_concentration_wet = open(particle_number_concentration_wet_file_path, "r+")
particle_number_concentration_wet = particle_number_concentration_wet.readlines()
particle_number_concentration_wet= pd.DataFrame([particle_number_concentration_wet[i].split(",") for i in range(2, len(particle_number_concentration_wet))])
particle_number_concentration_wet.to_csv('particle_number_concentration_wet_number.csv')

### Get information of tracked components
# get tracked components' names
model_variable_path = r"Model variable.txt"
model_variable = open(model_variable_path, "r+")
model_variable = model_variable.readlines()
model_variable = np.array([model_variable[i].split("=") for i in range(1, len(model_variable))])
tracked_comp = model_variable[model_variable[:,0] == "tracked_comp "]
tracked_comp = tracked_comp[0][1].split(", ")
tracked_comp[0] = tracked_comp[0][1:]
tracked_comp[-1] = tracked_comp[-1][:-1]
# get tracked components' indices
tracked_comp_index = [components.index(i) for i in tracked_comp if i in components]
# get tracked components' molecular weights
tracked_comp_weight = [molecular_weight[i] for i in tracked_comp_index]

for i in range(len(tracked_comp)):
    rate_of_change_file_path = tracked_comp[i] + "_rate_of_change"
    rate_of_change = open(rate_of_change_file_path, "r+")
    rate_of_change = rate_of_change.readlines()
    rate_of_change = np.array([rate_of_change[i].split(",") for i in range(1, len(rate_of_change))])
    rate_of_change = rate_of_change.astype(float)
    for j in range(len(rate_of_change[0])):
        rate_of_change[0][j] = rate_of_change[0][j] + 1

    # number concentration, unit: molecules/cc.s (air)
    rate_of_change_number = rate_of_change
    pd.DataFrame(rate_of_change_number).to_csv(tracked_comp[i]+"_rate_of_change_(molecules.(cc.s (air)).\u207B\u00b9).csv")

    # ppb
    rate_of_change_ppb = rate_of_change
    rate_of_change_ppb[1:] = (rate_of_change_ppb[1:].transpose() / factor[1:]).transpose()
    pd.DataFrame(rate_of_change_ppb).to_csv(tracked_comp[i]+"_rate_of_change (ppb.s\u207B\u00b9).csv")

    # mass concentration, unit: ug/m3
    temperature_new = np.tile(temperature[0], (len(rate_of_change) - 1, len(rate_of_change[0])))
    rate_of_change_mass = rate_of_change_ppb
    rate_of_change_mass[1:] = rate_of_change_mass[1:] * 12.187 / temperature_new * tracked_comp_weight[0]
    pd.DataFrame(rate_of_change_mass).to_csv(tracked_comp[i]+"_rate_of_change (\u03BCg.m\u207B\u00b3.s\u207B\u00b9).csv")

##### total_concentration_of_injected_components (mass concentration, unit: ug/m3)
total_concentration_of_injected_components_file_path = r"total_concentration_of_injected_components"
total_concentration_of_injected_components = open(total_concentration_of_injected_components_file_path, "r+")
total_concentration_of_injected_components = total_concentration_of_injected_components.readlines()
total_concentration_of_injected_components = np.array([total_concentration_of_injected_components[i].split(",") for i in range(1,len(total_concentration_of_injected_components))])
pd.DataFrame(total_concentration_of_injected_components).to_csv("total_concentration_of_injected_components (\u03BCg.m\u207B\u00b3).csv")


# Combine species information
species_information = pd.DataFrame({'Species':species, 
                                    'Molecular weight (g/mol)': molecular_weight, 
                                    'O:C ratio': OCratio,
                                    'Saturation vapor pressure at 298.15K':saturation_vapor_pressure})
organic_alkoxy_radical = [True if species_information.index[i] in organic_alkoxy_radical_index else False for i in range(len(species_information))]
species_information['Alkoxy Radicals'] = organic_alkoxy_radical
organic_peroxy_radical = [True if species_information.index[i] in organic_peroxy_radical_index else False for i in range(len(species_information))]
species_information['Peroxy Radicals'] = organic_peroxy_radical
species_information.to_csv('species_information.csv')

# Total particulate mass per species per size bins
perSpecies_perSizebin = pd.concat([pd.DataFrame(species), 
                                  pd.DataFrame(total_particulate_mass_species_sizebin)], 
                                  axis=1)
perSpecies_perSizebin.to_csv('perSpecies_perSizebin.csv')
                                 
# Total particulate mass per species per time
perSpecies_perTime = pd.concat([pd.DataFrame(species), 
                                  pd.DataFrame(total)], 
                                  axis=1)
perSpecies_perTime.to_csv('perSpecies_perTime.csv')

#size distribution SOA mass
particulate_phase_mass_df = pd.DataFrame(particulate_phase_mass)
particulate_phase_mass_df['Species'] = particulate_phase[0:particulate_phase.shape[0],0]
particulate_phase_mass_df['Size'] = particulate_phase[0:particulate_phase.shape[0],1]
#GroupBySize = particulate_phase_mass_df.groupby('Size')
#SizeMass = GroupBySize.sum()
#pd.DataFrame(SizeMass).to_csv('SOAsizemass.csv')
df = particulate_phase_mass_df[particulate_phase_mass_df.Species != ' core']
df = df[df.Species != ' H2O']
df = df[df.Species != ' AMMSUL']
pd.DataFrame(df).to_csv('SOAsizemass.csv')

#Other files
pd.DataFrame(gas_phase).to_csv('gas_phase.csv')
pd.DataFrame(concentration).to_csv("concentration_raw.csv")
