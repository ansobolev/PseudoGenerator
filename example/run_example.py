#!/usr/bin/env python

import numpy as np
import settings
from pseudogen.minimize_delta import minimize_delta
from pseudogen.calc_delta import read_ref_data


x0 = np.array([2.6, 2.8])
const_radii = ([2.3, 2.3, 0.7],)

print minimize_delta(settings, x0, const_radii)