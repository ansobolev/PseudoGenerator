#!/usr/bin/env python

import os
from log import get_logger, interlog
import uuid
import numpy as np
from generate import generate_pseudo, test_pseudo
from calc_delta import DeltaCalculation


def find_pseudo(settings):
    cwd = os.getcwd()
    element = settings.calc["element"]
    fdf_file = os.path.join(cwd, "siesta.fdf")
    # get uuid
    calc_uuid = uuid.uuid4().hex[:8]
    # make directories for the calculation
    dirname = os.path.join(element, calc_uuid)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    # logging
    logger = get_logger('find_pseudo', element)
    delta_calc = DeltaCalculation(settings, calc_uuid, logger)

    logger.info("Pseudo radii: {}".format(settings.radii))
    pseudo_file, err_pseudo = generate_pseudo(settings.calc, settings.electrons, settings.radii)
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
