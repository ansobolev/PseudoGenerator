#!/usr/bin/env python

""" Setuptools script for pseudogenerator
"""

from setuptools import setup, find_packages
from os.path import join, dirname

from pseudogen import __version__

setup(
    name='pseudogen',
    description='Siesta Pseudopotential Generator',
    version=__version__,
    packages=find_packages(),
    long_description=open(join(dirname(__file__), 'README.md')).read(),
)



