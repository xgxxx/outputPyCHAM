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

##### Size bin diameter (unit: um & mm)
diameter = radius*2
diameter_um = diameter
pd.DataFrame(diameter_um).to_csv('Diameter (\u03BCm).csv')

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

# nonSOA components, a txt file "nonSOA.txt" is given by user in which there is a line containing the names of all nonSOA components (e.g. core, H2O, AMMSUL)
nonSOA_path = r"nonSOA.txt"
nonSOA = open(nonSOA_path, "r+")
nonSOA = nonSOA.readlines()
nonSOA = nonSOA[0].split(", ")
components_number_SOA = components_number-len(nonSOA)

#Extract saturation vapor pressure at 298.15K of species
saturation_vapor_pressure = information[13].split(",")
saturation_vapor_pressure = saturation_vapor_pressure[1:len(saturation_vapor_pressure)]
saturation_vapor_pressure[0] = saturation_vapor_pressure[0][1:len(saturation_vapor_pressure[0])]
saturation_vapor_pressure[-1] = saturation_vapor_pressure[-1][0:len(saturation_vapor_pressure[-1])-2]
# Convert to int
saturation_vapor_pressure = [float(i) for i in saturation_vapor_pressure]

# Read chamber environment temperature
environment = open(chamber_environment_file_path, "r+")
environment = environment.readlines()
environment = environment[1:len(environment)]
temperature = [environment[i].split(",")[0] for i in range(len(environment))]
temperature = np.array(temperature).transpose().astype(float)
rh = [environment[i].split(",")[2] for i in range(len(environment))]
rh = np.array(rh).transpose().astype(float)

# Extract gas phase concentration (ppb)
gas_phase = concentration[concentration[:, 1] == "g"]
gas_phase[0, 0] = gas_phase[0, 0][2:]

# species corresponding to molecular weight
species = gas_phase[0:(len(gas_phase)), 0].tolist()

# Gas phase concentration for all components (ppb)
gas_phase_ppb_all = gas_phase
col_name = np.array([['Species', 'Phase'] + [str(i) + ' minute' for i in range(len(time))]])
gas_phase_ppb_all = np.vstack((col_name, gas_phase_ppb_all))
pd.DataFrame(gas_phase_ppb_all).to_csv("gas_phase_all (ppb).csv", index=False, header=False)

# Gas phase concentration for SOA components (ppb)
gas_phase_ppb_SOA = gas_phase_ppb_all
for i in nonSOA:
    gas_phase_ppb_SOA = gas_phase_ppb_SOA[gas_phase_ppb_SOA[:, 0] != ' ' + i]
pd.DataFrame(gas_phase_ppb_SOA).to_csv("gas_phase_SOA (ppb).csv", index=False, header=False)

# Gas phase concentration for nonSOA components (ppb)
gas_phase_ppb_nonSOA = np.array([i for i in gas_phase_ppb_all if i[0][1:] in nonSOA])
col_name = np.array([['Species', 'Phase'] + [str(i) + ' minute' for i in range(len(time))]])
gas_phase_ppb_nonSOA = np.vstack((col_name, gas_phase_ppb_nonSOA))
pd.DataFrame(gas_phase_ppb_nonSOA).to_csv("gas_phase_nonSOA (ppb).csv", index=False, header=False)

# repeat molecular weight array by no. of times of size bins
gas_phase = gas_phase[:, 2:]
gas_phase = gas_phase.astype(np.float)
gas_phase_mass = gas_phase * 12.187 / temperature
gas_phase_mass = (gas_phase_mass.transpose() * np.array(molecular_weight).transpose()).transpose()

# Gas phase concentration (ug/m3) for all components
gas_phase_mass_all = gas_phase_ppb_all[1:]
for i in range(len(gas_phase_mass_all)):
    gas_phase_mass_all[i] = np.concatenate((gas_phase_mass_all[i][:2], gas_phase_mass[i]))
col_name = np.array([['Species', 'Phase'] + [str(i) + ' minute' for i in range(len(time))]])
gas_phase_mass_all = np.vstack((col_name, gas_phase_mass_all))
pd.DataFrame(gas_phase_mass_all).to_csv("gas_phase_all (\u03BCg.m\u207B\u00b3).csv", index=False, header=False)

# Gas phase concentration (ug/m3) for SOA components
gas_phase_mass_SOA = gas_phase_mass_all
for i in nonSOA:
    gas_phase_mass_SOA = gas_phase_mass_SOA[gas_phase_mass_SOA[:, 0] != ' ' + i]
pd.DataFrame(gas_phase_mass_SOA).to_csv("gas_phase_SOA (\u03BCg.m\u207B\u00b3).csv", index=False, header=False)

# Gas phase concentration (ug/m3) for nonSOA components
gas_phase_mass_nonSOA = np.array([i for i in gas_phase_mass_all if i[0][1:] in nonSOA])
col_name = np.array([['Species', 'Phase'] + [str(i) + ' minute' for i in range(len(time))]])
gas_phase_mass_nonSOA = np.vstack((col_name, gas_phase_mass_nonSOA))
pd.DataFrame(gas_phase_mass_nonSOA).to_csv("gas_phase_nonSOA (\u03BCg.m\u207B\u00b3).csv", index=False, header=False)

# Gas phase concentration (molecules/cm3)
gas_phase_number = gas_phase * factor
# Gas phase concentration (molecules/cm3) for all components
gas_phase_number_all = gas_phase_mass_all[1:]
for i in range(len(gas_phase_number_all)):
    gas_phase_number_all[i] = np.concatenate((gas_phase_number_all[i][:2], gas_phase_number[i]))
col_name = np.array([['Species', 'Phase'] + [str(i) + ' minute' for i in range(len(time))]])
gas_phase_number_all = np.vstack((col_name, gas_phase_number_all))
pd.DataFrame(gas_phase_number_all).to_csv("gas_phase_all (molecules.(cc (air))\u207B\u00b9).csv", index=False, header=False)

# Gas phase concentration (molecules/cm3) for SOA components
gas_phase_number_SOA = gas_phase_number_all
for i in nonSOA:
    gas_phase_number_SOA = gas_phase_number_SOA[gas_phase_number_SOA[:, 0] != ' ' + i]
pd.DataFrame(gas_phase_number_SOA).to_csv("gas_phase_SOA (molecules.(cc (air))\u207B\u00b9).csv", index=False, header=False)

# Gas phase concentration (molecules/cm3) for nonSOA components
gas_phase_number_nonSOA = np.array([i for i in gas_phase_number_all if i[0][1:] in nonSOA])
col_name = np.array([['Species', 'Phase'] + [str(i) + ' minute' for i in range(len(time))]])
gas_phase_number_nonSOA = np.vstack((col_name, gas_phase_number_nonSOA))
pd.DataFrame(gas_phase_number_nonSOA).to_csv("gas_phase_nonSOA (molecules.(cc (air))\u207B\u00b9).csv", index=False, header=False)


# Extract particulate phase concentration (molecules/cm3)
particulate_phase = concentration[concentration[:, 1] != "g"]
particulate_phase = particulate_phase[particulate_phase[:, 1] != "w"]
wall_phase = concentration[concentration[:, 1] == "w"]
if len(wall_phase) == 0:
    particulate_phase[-1][1] = particulate_phase[-1][1][:-1]

# Particulate phase concentration (molecules/cm3) for all components
particulate_phase_all = particulate_phase
col_name = np.array([['Species', 'Sizebin'] + [str(i) + ' minute' for i in range(len(time))]])
particulate_phase_all = np.vstack((col_name, particulate_phase_all))
pd.DataFrame(particulate_phase_all).to_csv("particulate_phase_all (molecules.(cc (air))\u207B\u00b9).csv", index=False, header=False)

# Particulate phase concentration (molecules/cm3) for SOA components
particulate_phase_SOA = particulate_phase_all
for i in nonSOA:
    particulate_phase_SOA = particulate_phase_SOA[particulate_phase_SOA[:, 0] != ' ' + i]
pd.DataFrame(particulate_phase_SOA).to_csv("particulate_phase_SOA (molecules.(cc (air))\u207B\u00b9).csv", index=False, header=False)

# Particulate phase concentration (molecules/cm3) for nonSOA components
particulate_phase_nonSOA = np.array([i for i in particulate_phase_all if i[0][1:] in nonSOA])
col_name = np.array([['Species', 'Sizebin'] + [str(i) + ' minute' for i in range(len(time))]])
particulate_phase_nonSOA = np.vstack((col_name, particulate_phase_nonSOA))
pd.DataFrame(particulate_phase_nonSOA).to_csv("particulate_phase_nonSOA (molecules.(cc (air))\u207B\u00b9).csv", index=False, header=False)

# Convert from molecules/cm3 to ppb
particulate_phase_ppb = particulate_phase[:, 2:particulate_phase.shape[1]]
particulate_phase_ppb = particulate_phase_ppb.astype(np.float)
particulate_phase_ppb = particulate_phase_ppb / factor

# Particulate phase concentration (ppb) for all components
particulate_phase_ppb_all = particulate_phase
for i in range(len(particulate_phase)):
    particulate_phase_ppb_all[i] = np.concatenate((particulate_phase[i][:2], particulate_phase_ppb[i]))
col_name = np.array([['Species', 'Sizebin'] + [str(i) + ' minute' for i in range(len(time))]])
particulate_phase_ppb_all = np.vstack((col_name, particulate_phase_ppb_all))
pd.DataFrame(particulate_phase_ppb_all).to_csv("particulate_phase_all (ppb).csv", index=False, header=False)

# Particulate phase concentration (ppb) for SOA components
particulate_phase_ppb_SOA = particulate_phase_ppb_all
for i in nonSOA:
    particulate_phase_ppb_SOA = particulate_phase_ppb_SOA[particulate_phase_ppb_SOA[:, 0] != ' ' + i]
pd.DataFrame(particulate_phase_ppb_SOA).to_csv("particulate_phase_SOA (ppb).csv", index=False, header=False)

# Particulate phase concentration (ppb) for nonSOA components
particulate_phase_ppb_nonSOA = np.array([i for i in particulate_phase_ppb_all if i[0][1:] in nonSOA])
col_name = np.array([['Species', 'Sizebin'] + [str(i) + ' minute' for i in range(len(time))]])
particulate_phase_ppb_nonSOA = np.vstack((col_name, particulate_phase_ppb_nonSOA))
pd.DataFrame(particulate_phase_ppb_nonSOA).to_csv("particulate_phase_nonSOA (ppb).csv", index=False, header=False)

# repeat molecular weight array by no. of times of size bins
temp = np.tile(molecular_weight, int(len(particulate_phase_ppb) / len(molecular_weight)))
# Convert from ppb to ug/m3
particulate_phase_mass = particulate_phase_ppb * 12.187 / temperature
particulate_phase_mass = (particulate_phase_mass.transpose() * temp.transpose()).transpose()

# Particulate phase concentration (ug/m3) for all components
particulate_phase_mass_all = particulate_phase
for i in range(len(particulate_phase)):
    particulate_phase_mass_all[i] = np.concatenate((particulate_phase[i][:2], particulate_phase_mass[i]))
col_name = np.array([['Species', 'Sizebin'] + [str(i) + ' minute' for i in range(len(time))]])
particulate_phase_mass_all = np.vstack((col_name, particulate_phase_mass_all))
pd.DataFrame(particulate_phase_mass_all).to_csv("particulate_phase_all (\u03BCg.m\u207B\u00b3).csv", index=False, header=False)

# Particulate phase concentration (ug/m3) for SOA components
particulate_phase_mass_SOA = particulate_phase_mass_all
for i in nonSOA:
    particulate_phase_mass_SOA = particulate_phase_mass_SOA[particulate_phase_mass_SOA[:, 0] != ' ' + i]
pd.DataFrame(particulate_phase_mass_SOA).to_csv("particulate_phase_SOA (\u03BCg.m\u207B\u00b3).csv", index=False, header=False)

# Particulate phase concentration (ug/m3) for nonSOA components
particulate_phase_mass_nonSOA = np.array([i for i in particulate_phase_mass_all if i[0][1:] in nonSOA])
col_name = np.array([['Species', 'Sizebin'] + [str(i) + ' minute' for i in range(len(time))]])
particulate_phase_mass_nonSOA = np.vstack((col_name, particulate_phase_mass_nonSOA))
pd.DataFrame(particulate_phase_mass_nonSOA).to_csv("particulate_phase_nonSOA (\u03BCg.m\u207B\u00b3).csv", index=False, header=False)

### total mass concentration (ug/m3) per species per size bin
total_particulate_mass_species_sizebin = np.sum(particulate_phase_mass, axis=1)
total_particulate_mass_species_sizebin = total_particulate_mass_species_sizebin.reshape(bin_number, components_number).transpose()
# total mass concentration (ug/m3) per species per size bin for all components
total_particulate_mass_species_sizebin_all = particulate_phase[:components_number]
total_particulate_mass_species_sizebin_all = total_particulate_mass_species_sizebin_all[:, :bin_number + 1]
for i in range(components_number):
    total_particulate_mass_species_sizebin_all[i] = np.concatenate(
        (particulate_phase[i][:1], total_particulate_mass_species_sizebin[i]))
col_name = np.array([['Species'] + ['p ' + str(i + 1) for i in range(bin_number)]])
total_particulate_mass_species_sizebin_all = np.vstack((col_name, total_particulate_mass_species_sizebin_all))
pd.DataFrame(total_particulate_mass_species_sizebin_all).to_csv('perSpecies_perSizebin_all (\u03BCg.m\u207B\u00b3).csv', index=False, header=False)

# total mass concentration (ug/m3) per species per size bin for SOA components
total_particulate_mass_species_sizebin_SOA = total_particulate_mass_species_sizebin_all
for i in nonSOA:
    total_particulate_mass_species_sizebin_SOA = total_particulate_mass_species_sizebin_SOA[
        total_particulate_mass_species_sizebin_SOA[:, 0] != ' ' + i]
pd.DataFrame(total_particulate_mass_species_sizebin_SOA).to_csv('perSpecies_perSizebin_SOA (\u03BCg.m\u207B\u00b3).csv', index=False, header=False)

# total mass concentration (ug/m3) per species per size bin for nonSOA components
total_particulate_mass_species_sizebin_nonSOA = np.array(
    [i for i in total_particulate_mass_species_sizebin_all if i[0][1:] in nonSOA])
col_name = np.array([['Species'] + ['p ' + str(i + 1) for i in range(bin_number)]])
total_particulate_mass_species_sizebin_nonSOA = np.vstack((col_name, total_particulate_mass_species_sizebin_nonSOA))
pd.DataFrame(total_particulate_mass_species_sizebin_nonSOA).to_csv(
    'perSpecies_perSizebin_nonSOA (\u03BCg.m\u207B\u00b3).csv', index=False, header=False)

### total mass concentration (ug/m3) per species
total_particulate_mass_species = np.sum(total_particulate_mass_species_sizebin, axis=1)
# total mass concentration (ug/m3) per species for all components
total_particulate_mass_species_all = particulate_phase[:components_number, :2]
for i in range(components_number):
    total_particulate_mass_species_all[i] = np.concatenate(
        (particulate_phase[i][:1], [total_particulate_mass_species[i]]))
col_name = np.array([['Species', ' ']])
total_particulate_mass_species_all = np.vstack((col_name, total_particulate_mass_species_all))
pd.DataFrame(total_particulate_mass_species_all).to_csv('perSpecies_all (\u03BCg.m\u207B\u00b3).csv', index=False, header=False)

# total mass concentration (ug/m3) per species for SOA components
total_particulate_mass_species_SOA = total_particulate_mass_species_all
for i in nonSOA:
    total_particulate_mass_species_SOA = total_particulate_mass_species_SOA[total_particulate_mass_species_SOA[:, 0] != ' ' + i]
pd.DataFrame(total_particulate_mass_species_SOA).to_csv('perSpecies_SOA (\u03BCg.m\u207B\u00b3).csv', index=False, header=False)

# total mass concentration (ug/m3) per species for nonSOA components
total_particulate_mass_species_nonSOA = np.array([i for i in total_particulate_mass_species_all if i[0][1:] in nonSOA])
col_name = np.array([['Species', ' ']])
total_particulate_mass_species_nonSOA = np.vstack((col_name, total_particulate_mass_species_nonSOA))
pd.DataFrame(total_particulate_mass_species_nonSOA).to_csv('perSpecies_nonSOA (\u03BCg.m\u207B\u00b3).csv', index=False, header=False)

### total mass concentration (ug/m3) per species per time
total_particulate_mass_species_time = np.array_split(particulate_phase_mass, bin_number)
# total mass concentration (ug/m3) per species per time for all components
total = np.zeros(shape=(components_number, len(time)))
for items in total_particulate_mass_species_time:
    total = total + items
total_particulate_mass_species_time_all = particulate_phase[:components_number, :-1]
for i in range(components_number):
    total_particulate_mass_species_time_all[i] = np.concatenate((particulate_phase[i][:1], total[i]))
col_name = np.array([['Species'] + [str(i) + ' minute' for i in range(len(time))]])
total_particulate_mass_species_time_all = np.vstack((col_name, total_particulate_mass_species_time_all))
pd.DataFrame(total_particulate_mass_species_time_all).to_csv('perSpecies_perTime_all (\u03BCg.m\u207B\u00b3).csv', index=False, header=False)

# total mass concentration (ug/m3) per species per time for SOA components
total_particulate_mass_species_time_SOA = total_particulate_mass_species_time_all
for i in nonSOA:
    total_particulate_mass_species_time_SOA = total_particulate_mass_species_time_SOA[
        total_particulate_mass_species_time_SOA[:, 0] != ' ' + i]
pd.DataFrame(total_particulate_mass_species_time_SOA).to_csv('perSpecies_perTime_SOA (\u03BCg.m\u207B\u00b3).csv', index=False, header=False)

# total mass concentration (ug/m3) per species per time for nonSOA components
total_particulate_mass_species_time_nonSOA = np.array(
    [i for i in total_particulate_mass_species_time_all if i[0][1:] in nonSOA])
col_name = np.array([['Species'] + [str(i) + ' minute' for i in range(len(time))]])
total_particulate_mass_species_time_nonSOA = np.vstack((col_name, total_particulate_mass_species_time_nonSOA))
pd.DataFrame(total_particulate_mass_species_time_nonSOA).to_csv('perSpecies_perTime_nonSOA (\u03BCg.m\u207B\u00b3).csv', index=False, header=False)

# total SOA mass concentration for all components in each size bin at each time
total_sizebin_time_data = np.empty([bin_number, len(time)])
for bins in range(bin_number):
    temp_sum = particulate_phase_mass[components_number_SOA * bins:components_number_SOA * (bins + 1)]
    temp_sum = np.sum(temp_sum, axis=0)
    total_sizebin_time_data[bins] = temp_sum
total_sizebin_time = particulate_phase_mass_all[:bin_number+1, 1:]
for i in range(len(total_sizebin_time_data)):
    total_sizebin_time[i+1] = np.concatenate((['p'+str(i+1)], total_sizebin_time_data[i]))
pd.DataFrame(total_sizebin_time).to_csv('perSizebin perTime (\u03BCg.m\u207B\u00b3).csv', index=False, header=False)

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
# the path (path_model_var) of txt file "Model variable.txt" is given by user
path_model_var = r'C:\Users\24979\PyCHAM\PyCHAM\output\Chemical kpp file\41bins\inputs'
os.chdir(path_model_var)

# get tracked components' names
model_variable_path = r"Model variable.txt"
model_variable = open(model_variable_path, "r+")
model_variable = model_variable.readlines()
model_variable = np.array([model_variable[i].split("=") for i in range(1, len(model_variable))])
tracked_comp = model_variable[model_variable[:,0] == "tracked_comp "]
tracked_comp = tracked_comp[0][1].split(", ")
tracked_comp[0] = tracked_comp[0][1:]
tracked_comp[-1] = tracked_comp[-1][:-1]

os.chdir(path)
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
    rate_of_change_number = rate_of_change_number.transpose()
    pd.DataFrame(rate_of_change_number).to_csv(tracked_comp[i]+" rate of change (molecules.(cc.s (air)).\u207B\u00b9).csv", index=False, header=False)

    # ppb
    rate_of_change_ppb = rate_of_change
    rate_of_change_ppb[1:] = (rate_of_change_ppb[1:].transpose() / factor[:-1]).transpose()
    rate_of_change_ppb = rate_of_change_ppb.transpose()
    pd.DataFrame(rate_of_change_ppb).to_csv(tracked_comp[i]+" rate of change (ppb.s\u207B\u00b9).csv", index=False, header=False)

##### total_concentration_of_injected_components (mass concentration, unit: ug/m3)
total_concentration_of_injected_components_file_path = r"total_concentration_of_injected_components"
total_concentration_of_injected_components = open(total_concentration_of_injected_components_file_path, "r+")
total_concentration_of_injected_components = total_concentration_of_injected_components.readlines()
total_concentration_of_injected_components = np.array([total_concentration_of_injected_components[i].split(",") for i in range(1,len(total_concentration_of_injected_components))])
pd.DataFrame(total_concentration_of_injected_components).to_csv("total_concentration_of_injected_components (\u03BCg.m\u207B\u00b3).csv")

# O/C ratio
OCratio = information[7][8:-2].split(", ")
for i in range(len(OCratio)):
    o = OCratio[i].count('O')
    c = OCratio[i].count('C') + OCratio[i].count('c')
    OCratio[i] = o / c if c != 0 else 'no C'
OCratio = OCratio + [' '] + [' ']


# H/C ratio
HCratio = information[7][8:-2].split(", ")
for i in range(len(HCratio)):
    my_smiles_string = HCratio[i][1:-1]
    my_mol = Chem.MolFromSmiles(my_smiles_string)
    try:
        Chem.AddHs(my_mol)
    except:
        HCratio[i] = 'error'
    else:
        my_mol_with_explicit_h = Chem.AddHs(my_mol)
        h = my_mol_with_explicit_h.GetNumAtoms() - my_mol_with_explicit_h.GetNumHeavyAtoms()
        c = my_smiles_string.count('C') + my_smiles_string.count('c')
        HCratio[i] = h / c if c != 0 else 'no C'
HCratio = HCratio + [' '] + [' ']

# Combine species information
species_information = pd.DataFrame({'Species':species, 
                                    'Molecular weight (g/mol)': molecular_weight, 
                                    'O:C ratio': OCratio,
                                    'H:C ratio': HCratio,
                                    'Saturation vapor pressure at 298.15K':saturation_vapor_pressure})
organic_alkoxy_radical = [True if species_information.index[i] in organic_alkoxy_radical_index else False for i in range(len(species_information))]
species_information['Alkoxy Radicals'] = organic_alkoxy_radical
organic_peroxy_radical = [True if species_information.index[i] in organic_peroxy_radical_index else False for i in range(len(species_information))]
species_information['Peroxy Radicals'] = organic_peroxy_radical
species_information.to_csv('species_information.csv')

#Other files
pd.DataFrame(gas_phase).to_csv('gas_phase.csv')
pd.DataFrame(concentration).to_csv("concentration_raw.csv")
