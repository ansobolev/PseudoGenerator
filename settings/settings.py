# environment variables
ATOM_PROGRAM = '/home/ssozykin/distr/siesta/atom-4.0.2/atm'
ATOM_UTILS_DIR ='/home/ssozykin/distr/siesta/atom-4.0.2/Tutorial/Utils/'

element = "C" 

# general calculation parameters
calc = {"element": element,
        "alat": 2.468,
        "vectors": [[0.866025, -0.50000, 0.000000],
                    [0.000000,  1.00000, 0.000000],
                    [0.000000,  0.00000, 3.5818476]],
        "nat": 4,
        "lattice": "FCC",
        "xc": "pb",
        "n_core": 1,
        "n_val": 2,
        "is_spin_pol": False,
        "core": True, 
        }

# number of points to calculate BM EOS (odd value)
volumes = 3

# pseudopotential parameters
electrons = [2, 2]
radii = [1.54, 1.54, 1.54, 1.54, -1.0]

# SIESTA calculation parameters
siesta_calc = {"element": element,
               "title": element + " SIESTA calc",
               "xc_f": "GGA",
               "xc": "PBE"
              }

# electronic configurations
configs = [[1.5, 2.5],
           [2, 2],
           [2.5, 1.5],
           [1.25, 2.75],
           [3,1]]

# Minimization options

eps = 0.05
method = 'CG'
tolerance = 1e-2
min_options = {'disp': True}