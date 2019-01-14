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
from settings import *
from generate import PGInputFile, PTInputFile
from siesta import read_fdf_file, prepare_siesta_calc, run_siesta_calc
from get_energies import read_energy
from find_pseudo import siesta_calc
from calc_delta import BM, read_ref_data, calcDelta_one


if __name__ == "__main__":
    # fdf file
    fdf_file = read_fdf_file()
    data_dir = os.path.join("Fe", "37518ef0")
    #data_dir = os.path.join(element, sys.argv[1])
    
    x, y = [], []
    wd = os.getcwd()
    os.chdir(data_dir)
    pseudo_file = glob.glob("*.psf")[0]
    for root, dirs, files in os.walk(os.getcwd()):
        if "check" in root: continue
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
    x = (np.array(x) ** 3) / nat
    y = np.array(y)
    p = np.polyfit(x, y, 2)
    # find equilibrium volume
    x_min = -p[1] / (2*p[0])
    alats = (np.linspace(0.94*x_min, 1.06*x_min, 7) * nat) ** (1./3.)
    np.savetxt("text.txt", np.vstack((alats, np.polyval(p,np.linspace(0.94*x_min, 1.06*x_min, 7)))).T)
    # get check directory
    if not os.path.exists("check"):
        os.makedirs("check")
    shutil.copy(pseudo_file, "check")
    os.chdir("check")
    x = []
    y = []
    for alat in alats:
      # prepare_siesta_calc(fdf_file, pseudo_file, alat, siesta_calc)
      # run_siesta_calc(alat, siesta_calc)
        print alat
        e = read_energy(alat)
        if e is not None:
            x_i, y_i = e
            x.append(x_i)
            y.append(y_i)
    x_p = (np.array(x) ** 3) / nat
    y_p = np.array(y)
    print x, y
    # print x_p, y_p
    vol, bulk_mod, bulk_deriv, res = BM(np.vstack((x_p, y_p)).T)
    np.savetxt("energies.txt", np.vstack((x_p, y_p)).T)
    
    # our_data = np.core.records.fromrecords([(siesta_calc['element'], vol, bulk_mod, bulk_deriv), ], names=('element', 'V0', 'B0', 'BP'))
    # delta, delta_rel, delta1 = calcDelta_one(our_data, ref_data_el, useasymm=False)
    # log["delta"] = delta
    # log["delta_rel"] = delta_rel
    # log["delta1"] = delta1
    # log_entry(log, calc["element"])


