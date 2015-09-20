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
from settings import element, nat
from delta.eosfit import BM
from delta.calcDelta import calcDelta_one, read_ref_data

def read_energy(alat):
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
            if "Final energy" in line:
                final = True
            if final and "Total" in line:
                energy = float(line.split()[3])
                final = False
                x = float(path.split("/")[-1])
                y = energy
                return np.array((x, y))
    return None


if __name__ == "__main__":
    # if len(sys.argv) < 2:
    #     raise SystemExit, "Usage: get_energies.py <data directory>"
    
    # data_dir = sys.argv[1]
    data_dir = "Al/4ce268f7"
    x = []
    y = []
    wd = os.getcwd()
    # os.chdir(os.path.join(data_dir, "check"))
    os.chdir(data_dir)
    for root, dirs, files in os.walk(os.getcwd()):
        for dir_i in dirs:
            try:
                alat = float(dir_i)
            except:
                continue
            energies = read_energy(alat)
            if energies is not None:
                x_i, y_i = energies
                x.append(x_i)
                y.append(y_i)
    os.chdir(wd)
    x = (np.array(x) ** 3) / nat
    y = np.array(y)
    p = np.polyfit(x, y, 2)
    x_min = -p[1] / (2*p[0])
    x_p = np.linspace(0.94*x_min, 1.06*x_min, 7)
    y_p = np.polyval(p, x_p)
#    vol, bulk_mod, bulk_deriv, res = BM(np.vstack((x_p, y_p)).T)
    vol, bulk_mod, bulk_deriv, res = BM(np.vstack((x, y)).T)
    ref_data = read_ref_data("delta/WIEN2k.txt")
    ref_data_el = ref_data[ref_data['element'] == element]
    our_data = np.core.records.fromrecords([(element, vol, bulk_mod, bulk_deriv), ], names=('element', 'V0', 'B0', 'BP'))

    print calcDelta_one(our_data, ref_data_el, useasymm=False)
    
    plt.plot(np.array(x), np.array(y), "o")
    plt.plot(x_p, y_p, "--")
    plt.show()

