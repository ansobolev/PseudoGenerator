# environment variables
ATOM_PROGRAM = '/home/physics/bin/atm'
ATOM_UTILS_DIR ='/home/physics/bin/pseudo'

element = "Al"

equil_volume = 16.4796

# general calculation parameters
calc = {"element": element,
        "lattice": "FCC",
        "xc": "pb",
        "n_core": 3,
        "n_val": 2,
        "is_spin_pol": False,
        "core": True, 
        }

# pseudopotential parameters
electrons = [2, 1]
radii = [2.4, 2.8, 2.3, 2.3, 0.7]

# SIESTA calculation parameters
siesta_calc = {"element": element,
               "title": element + " SIESTA calc",
               "xc_f": "GGA",
               "xc": "PBE"
              }

# electronic configurations
configs = [[1.5, 1.5],
           [1, 2],
           [0.5, 2.5],
           [0, 3]]

# number of atoms in cubic cell
_nat_cell = {"SC": 1,
             "BCC": 2,
             "FCC": 4}
nat = _nat_cell[calc["lattice"]]
