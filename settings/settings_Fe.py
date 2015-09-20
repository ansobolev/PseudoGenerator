# environment variables
ATOM_PROGRAM = '/home/physics/bin/atm'
ATOM_UTILS_DIR ='/home/physics/bin/pseudo'

element = "Fe"

equil_volume = 11.3436

# general calculation parameters
calc = {"element": element,
        "lattice": "BCC",
        "xc": "pb",
        "n_core": 5,
        "n_val": 2,
        "is_spin_pol": True,
        "core": True, 
        }

# pseudopotential parameters
electrons = [(2,0), (6,0)]
radii = [2., 2.25, 2., 2., 0.7]

# SIESTA calculation parameters
siesta_calc = {"element": element,
               "title": element + " SIESTA calc",
               "xc_f": "GGA",
               "xc": "PBE"
              }

# electronic configurations
configs = [[(1,0),(7,0)], 
           [(2,0),(6,0)],
           [(1.5,0),(6.5,0)],
           [(0,0),(8,0)],
           [(1,0),(6,0),(1,0)],
           [(1,0),(5,0),(2,0)]]

# number of atoms in cubic cell
_nat_cell = {"SC": 1,
             "BCC": 2,
             "FCC": 4}
nat = _nat_cell[calc["lattice"]]