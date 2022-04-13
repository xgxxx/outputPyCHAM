''' Module to plot contour for SOA formed over time and size bin  '''

import numpy as np
import matplotlib.pyplot as plt


def contour_SOA(time, bin_number, components_number_SOA, particulate_phase_mass_SOA, resolution):
    # time: given in post processing file, which is extracted from PyCHAM output time file
    # bin_number: number of size bins, given in post processing file
    # components_number_SOA: number of SOA components, given in post processing file
    # particulate_phase_mass_SOA: contains partitulate phase mass concentration\
    # for each SOA component in each size bin and each time, given in post processing file
    # resolution: determine number of lines in contour plot, given by user

    ###x: time (unit:min)
    x = time

    ###y: size bin no.
    y = np.arange(1, bin_number+1, 1)

    ###z: total mass concentration (ug/m3) per size bin
    z = np.zeros((bin_number, len(time)))
    for i in range(bin_number):
        # calculate total SOA mass concentration in each size bin at each time
        t = particulate_phase_mass_SOA[components_number_SOA * i:components_number_SOA * (i+1)]
        z[i] = np.sum(t, axis=0)

    ###plot contour
    plt.figure(1)
    plt.contourf(x, y, z, resolution, cmap='RdBu')
    plt.colorbar(label='unit: \u03BCg/m\u00b3')
    plt.title('Contour for SOA (\u03BCg/m\u00b3) formed over time and size bin')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Size bin number')
    plt.savefig('Contour for SOA formed over time and size bin')
