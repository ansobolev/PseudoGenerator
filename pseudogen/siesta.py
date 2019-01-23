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

def read_fdf_file(file_name):
    with open(file_name, "r") as f:
        file_text = f.read()
    return file_text

def write_fdf_file(file_name, file_text, calc):
    # write vectors in a cool way
    calc["vectors"] = "\n".join(["\t".join([str(v) for v in calc["vectors"]])])
    with open(file_name, "w") as f:
        f.write(file_text.format(**calc))

class SiestaCalculation(object):

    def __init__(self, settings, pseudo_file, fdf_file="siesta.fdf"):
        self.calc = settings.calc
        self.siesta_calc = settings.siesta_calc
        self.element = self.calc["element"]
        
        self.pseudo_file = pseudo_file  
        self.fdf_file = read_fdf_file(fdf_file)

    def prepare(self, alat=None):
        if alat is not None:
            self.siesta_calc["alat"] = alat
        else:
            self.siesta_calc["alat"] = self.calc["alat"]
        self.siesta_calc["vectors"] = self.calc["vectors"]
        path = "%.4f" % (alat,)
        if not os.path.isdir(path):
            os.makedirs(path)
        fdf_file_name = os.path.join(path, self.element + ".fdf")
        # copy and rename psf file
        shutil.copy(self.pseudo_file, path)
        os.rename(os.path.join(path, self.pseudo_file), os.path.join(path, self.element + ".psf"))
        write_fdf_file(fdf_file_name, self.fdf_file, self.siesta_calc)

    def run(self):
        path = "%.4f" % (self.siesta_calc["alat"],)
        cwd = os.getcwd()
        os.system('cd ' + path + '; mpirun -np 4 ' + SIESTA_EXEC + ' < ' + self.element + ".fdf" + ' | tee ' + self.element + '.out')
        os.chdir(cwd)

    def results(self):
        out_file = self.element + '.out'
        # get files in path
        path = "%.4f" % (self.siesta_calc["alat"],)
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