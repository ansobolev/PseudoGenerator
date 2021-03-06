# Siesta Pseudopotential Generator

The generator of ''good'' pseudopotentials for SIESTA code. To make it work, several things have to be taken care of:

 1. ATOM program (from `Pseudo` directory in SIESTA) has to be installed in the system 
 2. `settings.py` file. It can be taken from `settings` folder and renamed accordingly; 
or it can be written from scratch using any file from this folder as an example.
 3.  files from **[DeltaCodesDFT]** project (actually, only `delta/WIEN2k.txt` is needed, to calculate delta factor against reference code). These files must be placed in the calculation folder in `delta` subdir.

Also, several environmental variables have to be defined in the system:
 
 * `ATOM_PROGRAM`, which defines the location of the ATOM `atm` executable;
 * `ATOM_UTILS_DIR`, which defines the location of the `pt.sh`, `pg.sh` and `ae.sh` scripts;
 * `SIESTA_EXEC`, which defines the location of `SIESTA` executable.

## Installation

First, clone the code from GitHub:

```
git clone https://github.com/ansobolev/PseudoGenerator.git
```

Then, install the package with `pip`:
```
pip install -e ./PseudoGenerator
```

## Usage 

Look at the `examples` folder for an overview of the usage.

[DeltaCodesDFT]: <http://molmod.ugent.be/deltacodesdft>
