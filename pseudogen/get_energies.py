#!/usr/bin/env python

"""
Plots total energies vs volume for directory given as a command-line argument
"""

__author__ = "Andrey Sobolev"
__year__ = "2015"

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from calc_delta import BM, read_ref_data, calcDelta


def read_energy(element, alat):
    out_file = element + '.out'
    # get files in path
    path = "%.4f" % (alat,)
    files = os.listdir(path)
    # check if the calc succeeded 
    if out_file in files:
        with open(os.path.join(path, out_file), "r") as f:
            lines = f.readlines()
        final = False
        for line in lines:
            if "Cell volume" in line:
                vol = float(line.split()[-1])
            if "Final energy" in line:
                final = True
            if final and "Total" in line:
                energy = float(line.split()[3])
                final = False
                return np.array((vol, energy))
    return None


def get_energies(settings, data_dir):
    x = []
    y = []
    wd = os.getcwd()
    element = settings.calc["element"]
    # os.chdir(os.path.join(data_dir, "check"))
    os.chdir(data_dir)
    for _, dirs, _ in os.walk(os.getcwd()):
        for dir_i in dirs:
            try:
                alat = float(dir_i)
            except:
                continue
            energies = read_energy(element, alat)
            if energies is not None:
                x_i, y_i = energies
                x.append(x_i)
                y.append(y_i)
    os.chdir(wd)
    x = (np.array(x) ** 3) / settings.nat
    y = np.array(y)
    p = np.polyfit(x, y, 2)
    x_min = -p[1] / (2*p[0])
    x_p = np.linspace(0.94*x_min, 1.06*x_min, 7)
    y_p = np.polyval(p, x_p)
#    vol, bulk_mod, bulk_deriv, res = BM(np.vstack((x_p, y_p)).T)
    vol, bulk_mod, bulk_deriv, _ = BM(np.vstack((x, y)).T)
    ref_data = read_ref_data("delta/WIEN2k.txt")
    ref_data_el = ref_data[ref_data['element'] == element]
    our_data = np.core.records.fromrecords([(element, vol, bulk_mod, bulk_deriv), ], names=('element', 'V0', 'B0', 'BP'))

    print calcDelta(our_data, ref_data_el, useasymm=False)
    
    plt.plot(np.array(x), np.array(y), "o")
    plt.plot(x_p, y_p, "--")
    plt.show()

