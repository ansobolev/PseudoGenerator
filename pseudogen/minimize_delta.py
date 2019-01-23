#!/usr/bin/env python

"""
minimize_delta.py minimizes delta-factor by conjugate gradients method
"""

import os
import uuid
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from generate import generate_pseudo, test_pseudo
from siesta import SiestaCalculation
from get_energies import read_energy
from calc_delta import BM, read_ref_data, calcDelta

def log_entry(log, element="Fe"):
    entry = "   " + "-" * 10 + "   "
    entry += """
Dir = {uuid}    Radii = {radii}
Pseudo error (ground state) = {err_pseudo:.4} Ry
                                  max      mean
Pseudo error (test configs) =    {err_max:6.4}  {err_mean:6.4} Ry
Equilibrium volume per atom =    {min_p:6.6} A^3
                                  delta  rel_delta
Delta factor                =    {delta:6.4}  {delta_rel:6.4} meV/atom
             """
    with open(element + "/log.dat", "a") as f:
        f.write(entry.format(**log))


def minimize_delta(settings, x0, const_radii):
    cwd = os.getcwd()
    fdf_file = os.path.join(cwd, "siesta.fdf")
    element = settings.calc["element"]
    # get reference data
    ref_data = read_ref_data("delta/WIEN2k.txt")
    ref_data_el = ref_data[ref_data['element'] == element]

    def fun(args, consts):
        log = {}
        # get uuid
        calc_uuid = uuid.uuid4().hex[:8]
        # logging
        log["uuid"] = calc_uuid
        # make directories for the calculation
        dirname = os.path.join(element, calc_uuid)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        os.chdir(dirname)
        # generate radii iterable 
        radii = list(args) + list(consts)
        log["radii"] = radii
        pseudo_file, log["err_pseudo"] = generate_pseudo(settings.calc, settings.electrons, radii)
        log["err_mean"], log["err_max"] = test_pseudo(settings.calc, settings.configs)
        volumes = np.linspace(0.98, 1.02, 3)
        alats = (volumes * settings.equil_volume * settings.nat) ** (1./3.)
        x, y = [], []
        siesta_calc = SiestaCalculation(settings, pseudo_file, fdf_file=fdf_file)
        for alat in alats:
            siesta_calc.prepare(alat)
            siesta_calc.run()
            e = siesta_calc.results()
            if e is not None:
                x.append(float(e[0]))
                y.append(e[1])
        os.chdir(cwd)
        # making x and y arrays
        x = np.array(x) / settings.nat
        y = np.array(y) / settings.nat
        p = np.polyfit(x, y, 2)
        min_p = -p[1] / (2*p[0])
        log["min_p"] = min_p
        x_p = np.linspace(0.94*min_p, 1.06*min_p, 7)
        y_p = np.polyval(p, x_p)
        vol, bulk_mod, bulk_deriv, _ = BM(np.vstack((x_p, y_p)).T)
        our_data = np.core.records.fromrecords([(element, vol, bulk_mod, bulk_deriv), ], names=('element', 'V0', 'B0', 'BP'))
        delta, delta_rel, delta1 = calcDelta(our_data, ref_data_el, useasymm=False)
        log["delta"] = delta
        log["delta_rel"] = delta_rel
        log["delta1"] = delta1
        log_entry(log, element)
        return delta

    return minimize(fun, x0, args=const_radii, options={"eps":0.1})
