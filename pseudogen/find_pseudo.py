#!/usr/bin/env python

import os
import uuid
import numpy as np
import matplotlib.pyplot as plt
from generate import generate_pseudo, test_pseudo
from siesta import read_fdf_file, prepare_siesta_calc, run_siesta_calc
from get_energies import read_energy
from calc_delta import BM, read_ref_data, calcDelta


def log_entry(log, element):
    entry = "   " + "-" * 10 + "   "
    entry += """
Dir = {uuid}    Radii = {radii}
Pseudo error (ground state) = {err_pseudo:.4} Ry
                                  max      mean
Pseudo error (test configs) =    {err_max:6.4}  {err_mean:6.4} Ry
Equilibrium volume per atom =    {min_p:6.6} A^3
                                  delta, meV/atom  rel_delta, %
Delta factor                =    {delta:6.4}       {delta_rel:6.4}
             """
    with open(element + "/log.dat", "a") as f:
        f.write(entry.format(**log))


def get_volumes(n_vol, calc, alat=None):
    abc = np.array(calc["vectors"])
    vol = np.linalg.det(abc) * (calc["alat"] ** 3)
    k = (n_vol - 1) / 2
    return np.linspace(1-0.02*k, 1+0.02*k, n_vol) * vol

def get_alats(volumes, calc):
    vol_abc = np.linalg.det(np.array(calc["vectors"]))
    return (volumes / vol_abc) ** (1./3)

def find_pseudo(settings):
    cwd = os.getcwd()
    fdf_file = read_fdf_file()
    log = {}
    element = settings.calc["element"]
    # get uuid
    calc_uuid = uuid.uuid4().hex[:8]
    # logging
    log["uuid"] = calc_uuid
    # get reference data
    ref_data = read_ref_data("delta/WIEN2k.txt")
    ref_data_el = ref_data[ref_data['element'] == element]
    # make directories for the calculation
    dirname = os.path.join(element, calc_uuid)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    os.chdir(dirname)
    log["radii"] = settings.radii
    pseudo_file, log["err_pseudo"] = generate_pseudo(settings.calc, settings.electrons, settings.radii)
    log["err_mean"], log["err_max"] = test_pseudo(settings.calc, settings.configs)
    
    # possible alats (several points near Wien2K equilibrium)
    volumes = get_volumes(settings.volumes, settings.calc)
    alats = get_alats(volumes, settings.calc)
    x, y = [], []
    for alat in alats:
        prepare_siesta_calc(fdf_file, pseudo_file, alat, settings.siesta_calc)
        run_siesta_calc(alat, settings.siesta_calc)
        e = read_energy(element, alat)
        if e is not None:
            x.append(float(e[0]))
            y.append(e[1])
    os.chdir(cwd)
    # making x and y arrays
    x = np.array(x) / settings.calc["nat"]
    y = np.array(y) / settings.calc["nat"]
    p = np.polyfit(x, y, 2)
    min_p = -p[1] / (2*p[0])
    log["min_p"] = min_p
    x_p = np.linspace(0.94*min_p, 1.06*min_p, 7)
    y_p = np.polyval(p, x_p)
    vol, bulk_mod, bulk_deriv, _ = BM(np.vstack((x_p, y_p)).T)
    our_data = np.core.records.fromrecords([(element, vol, bulk_mod, bulk_deriv), ], names=('element', 'V0', 'B0', 'BP'))
    print our_data
    delta, delta_rel, delta1 = calcDelta(our_data, ref_data_el, useasymm=False)
    log["delta"] = delta
    log["delta_rel"] = delta_rel
    log["delta1"] = delta1
    log_entry(log, element)
