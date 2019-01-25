#!/usr/bin/env python

""" check_pseudo.py calculates energy for 7 alat points near SIESTA equilibrium to fine tune the delta-factor. 
"""

import os
import sys
import uuid
import glob
import numpy as np
import shutil
import matplotlib.pyplot as plt
from generate import PGInputFile, PTInputFile
from get_energies import read_energy
from calc_delta import BM, read_ref_data, calcDelta, get_alats, get_volumes


def check_pseudo(settings, data_dir):
    """ Checks pseudopotential for delta factor calculation
    
    Arguments:
        settings {[type]} -- [description]
        data_dir {[type]} -- [description]
    """
    cwd = os.getcwd()
    element = settings.calc["element"]
    x, y = [], []
    os.chdir(data_dir)
    pseudo_file = glob.glob("*.psf")[0]
    
    for root, dirs, _ in os.walk(os.getcwd()):
        if "check" in root: continue
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
    x = np.array(x) / settings.calc["nat"]
    y = np.array(y) / settings.calc["nat"]
    p = np.polyfit(x, y, 2)

    # make 7 points out of existing data
    if len(x) == 7:
        x_p = x
        y_p = y
    else:
        x_p = get_volumes(7, settings.calc) / settings.calc["nat"]
        y_p = np.poly1d(p)(x_p)

    # get check directory
    if not os.path.exists("check"):
        os.makedirs("check")
    shutil.copy(pseudo_file, "check")
    os.chdir("check")

    # write original data
    np.savetxt("energies_original.txt", np.vstack((x, y)).T)

    vol, bulk_mod, bulk_deriv, _ = BM(np.vstack((x_p, y_p)).T)
    np.savetxt("energies_BM.txt", np.vstack((x_p, y_p)).T)
    
    our_data = np.core.records.fromrecords([(element, vol, bulk_mod, bulk_deriv), ], names=('element', 'V0', 'B0', 'BP'))
    ref_data = read_ref_data(os.path.join(cwd, "delta", "WIEN2k.txt"))
    ref_data_el = ref_data[ref_data['element'] == element]
    delta, delta_rel, _ = calcDelta(our_data, ref_data_el, useasymm=False)
    with open("BP.dat", "w") as f:
        f.write("Our data: {}\n".format(our_data))
        f.write("Reference data: {}\n".format(ref_data_el))
        f.write("Delta factor: {} {}\n".format(delta, delta_rel))
