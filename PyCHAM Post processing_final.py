#####################################
## BN ##
import numpy as np
import pandas as pd 

# Set working directory - input folder of interest
import os
path = r'C:\Users\janet\Documents\nBox\0. SOAFP-PyCHAM output\PR3_simulation\Simulation results\20Nov_1700_1759_89VOCs_20bins'
os.chdir(path)

##### Size bin radius
radius_file_path = r"size_bin_radius"
radius = open(radius_file_path, "r+")
radius = radius.readlines()
radius_mm = pd.DataFrame([radius[i].split(",") for i in range(2, len(radius))])
radius_mm.to_csv('Radius (mm).csv')

##### Concentration information
# Input filepath here #
in_file_path = r"concentrations_all_components_all_times_gas_particle_wall"
information_file_path = r"model_and_component_constants"
chamber_environemnt_file_path = r"chamber_environmental_conditions"

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
# O/C ratio
OCratio = information[14].split(",")
# Fix the first
OCratio = OCratio[1:len(OCratio)] 
OCratio[0] = OCratio[0][2:len(OCratio[0])]
OCratio[-1] = OCratio[-1][0:len(OCratio[-1])-3]
# Convert to float
OCratio = [float(i) for i in OCratio]

# Read chamber environment temperature
environment = open(chamber_environemnt_file_path, "r+")
environment = environment.readlines()
environment = environment[1:len(environment)]
temperature = [environment[i].split(",")[0] for i in range(len(environment))]
temperature = np.array(temperature).transpose().astype(float)
rh = [environment[i].split(",")[2] for i in range(len(environment))]
rh = np.array(rh).transpose().astype(float)

# Extract particulate phase concentration
particulate_phase = concentration[concentration[:,1] != "g"]
# Convert from molecules/cm3 to ppb
particulate_phase_ppb = particulate_phase[:,2:particulate_phase.shape[1]]
particulate_phase_ppb = particulate_phase_ppb.astype(np.float)
particulate_phase_ppb = particulate_phase_ppb / factor
# repeat molecular weight array by no. of times of size bins
temp = np.tile(molecular_weight,int(len(particulate_phase_ppb)/len(molecular_weight)))
# Convert from ppb to ug/m3
particulate_phase_mass = particulate_phase_ppb*12.187/temperature
particulate_phase_mass = (particulate_phase_mass.transpose()*temp.transpose()).transpose()

# total mass concentration (ug/m3) per species per size bin
total_particulate_mass_species_sizebin = np.sum(particulate_phase_mass, axis=1)
total_particulate_mass_species_sizebin = total_particulate_mass_species_sizebin.reshape(int(len(particulate_phase_ppb)/len(molecular_weight)),len(molecular_weight)).transpose()
# total mass concentration (ug/m3) per species
total_particulate_mass_species = np.sum(total_particulate_mass_species_sizebin, axis = 1)
# total mass concentration (ug/m3) per species per time
total_particulate_mass_species_time = np.array_split(particulate_phase_mass, int(len(particulate_phase_ppb)/len(molecular_weight)))
total = np.zeros(shape=(len(molecular_weight),len(environment)))
for items in total_particulate_mass_species_time:
    total = total + items

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
