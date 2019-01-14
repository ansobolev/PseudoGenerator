#!/usr/bin/env python

"""
The script takes example fdf file and fills it with values, then makes subdirectories under 
the {element} folder corresponding to several per-atom volumes and runs SIESTA calculations
in each of the subdirectories. Pseudopotential must be provided as command-line argument.

The file is licensed under MIT License.
"""

__author__ = "Andrey Sobolev"
__year__ = 2015

import os
import sys
import shutil
import subprocess
import numpy as np

SIESTA_EXEC = os.environ.get('SIESTA_EXEC', '/home/andrey/bin/siesta')

def read_fdf_file(file_name="siesta.fdf"):
    with open(file_name, "r") as f:
        file_text = f.read()
    return file_text

def write_fdf_file(file_name, file_text, calc):
    with open(file_name, "w") as f:
        f.write(file_text.format(**calc))

def prepare_siesta_calc(fdf_file, pseudo_file, alat, calc):
    calc["alat"] = alat
    path = "%.4f" % (alat,)
    if not os.path.isdir(path):
        os.makedirs(path)
    fdf_file_name = os.path.join(path, calc["element"] + ".fdf")
    # copy and rename psf file
    shutil.copy(pseudo_file, path)
    os.rename(os.path.join(path, pseudo_file), os.path.join(path, calc["element"] + ".psf"))
    write_fdf_file(fdf_file_name, fdf_file, calc)

def run_siesta_calc(alat, calc):
    path = "%.4f" % (alat,)
    cwd = os.getcwd()
    os.system('cd ' + path + '; mpirun -np 4 ' + SIESTA_EXEC + ' < ' + calc["element"] + ".fdf" + ' | tee ' + calc["element"] + '.out')
    os.chdir(cwd)

