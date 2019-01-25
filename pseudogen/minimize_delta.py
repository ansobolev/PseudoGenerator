#!/usr/bin/env python

"""
minimize_delta.py minimizes delta-factor by conjugate gradients method
"""

import os
import uuid
import numpy as np
from scipy.optimize import minimize
from generate import generate_pseudo, test_pseudo
from calc_delta import DeltaCalculation
from log import get_logger, interlog

def minimize_delta(settings, x0, const_radii):
    cwd = os.getcwd()
    print cwd
    fdf_file = os.path.join(cwd, "siesta.fdf")
    element = settings.calc["element"]

    def fun(args, consts):
        # get uuid
        calc_uuid = uuid.uuid4().hex[:8]
        # make directories for the calculation
        dirname = os.path.join(element, calc_uuid)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        # logging
        logger = get_logger('find_pseudo', element)
        delta_calc = DeltaCalculation(settings, calc_uuid, logger)
        os.chdir(os.path.join(cwd, dirname))
        # generate radii iterable
        radii = list(args) + list(consts)
        logger.info("Pseudo radii: {}".format(radii))
        pseudo_file, err_pseudo = generate_pseudo(settings.calc, settings.electrons, radii)
        err_mean, err_max = test_pseudo(settings.calc, settings.configs)
        message = """
        Pseudo error (ground state) = {err_pseudo:.4} Ry
                                        max      mean
        Pseudo error (test configs) =    {err_max:6.4}  {err_mean:6.4} Ry"""
        logger.info(message.format(err_pseudo=err_pseudo, err_max=err_max, err_mean=err_mean))
        delta_calc.add_pseudo(pseudo_file)
        delta_calc.run_calcs(fdf_file)
        os.chdir(cwd)
        delta_calc.get_delta()
        interlog(logger)
        return delta_calc.delta

    eps = getattr(settings, 'eps', 0.1)
    return minimize(fun, x0, args=const_radii, options={"eps":eps})
