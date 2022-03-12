''' Module to plot contour for SOA formed over time and size bin  '''

import numpy as np
import matplotlib.pyplot as plt


def contour_SOA(time, radius, bin_number, components_number, particulate_phase_mass):
    ###x: time
    x = time

    ###y: size_bin_radius
    y = radius

    ###z: total mass concentration (ug/m3) per size bin
    z = np.zeros((bin_number, len(time)))
    for i in range(bin_number):
        t = particulate_phase_mass[components_number * i:components_number * (i+1)]
        z[i] = np.sum(t, axis=0)
    print(z)

    ###plot contour
    plt.figure(1)
    plt.contourf(x, y, z)
    plt.colorbar()
    plt.title("Contour for SOA (ug/m3) formed over time and size bin")
    plt.xlabel('time (min)')
    plt.ylabel('size bin radius (um)')
    plt.savefig("Contour for SOA formed over time and size bin")

