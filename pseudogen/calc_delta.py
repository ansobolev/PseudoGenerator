"""
 The following file is the set of routines taken and adapted from DeltaCodesDFT project,
 https://molmod.ugent.be/deltacodesdft
 Copyright (C) 2012 Kurt Lejaeghere <Kurt.Lejaeghere@UGent.be>, Center for
 Molecular Modeling (CMM), Ghent University, Ghent, Belgium
"""

import os
import numpy as np
from siesta import SiestaCalculation

def read_ref_data(file_name):
    """ Read reference data on V0, B0 and B1 from file
    
    Arguments:
        file_name {string} -- name of the file containing data
    """
    return np.genfromtxt(file_name, 
                         names=("element", "V0", "B0", "BP"), 
                         comments="#", 
                         dtype=None)



def BM(energies):

    fitdata = np.polyfit(energies[:,0]**(-2./3.), energies[:,1], 3, full=True)
    ssr = fitdata[1]
    sst = np.sum((energies[:,1] - np.average(energies[:,1]))**2.)
    residuals0 = ssr/sst
    deriv0 = np.poly1d(fitdata[0])
    deriv1 = np.polyder(deriv0, 1)
    deriv2 = np.polyder(deriv1, 1)
    deriv3 = np.polyder(deriv2, 1)

    volume0 = 0
    x = 0
    for x in np.roots(deriv1):
        if x > 0 and deriv2(x) > 0:
            volume0 = x**(-3./2.)
            break

    if volume0 == 0:
        print('Error: No minimum could be found')
        exit()
    
    derivV2 = 4./9. * x**5. * deriv2(x)
    derivV3 = (-20./9. * x**(13./2.) * deriv2(x) -
        8./27. * x**(15./2.) * deriv3(x))
    bulk_modulus0 = derivV2 / x**(3./2.)
    bulk_deriv0 = -1 - x**(-3./2.) * derivV3 / derivV2

    return volume0, bulk_modulus0, bulk_deriv0, residuals0

def calcDelta(data_f, data_w, useasymm):
    """
    Calculate the Delta using the data in data_f, data_w 
    """

    v0w = data_w['V0']
    b0w = data_w['B0'] * 10.**9. / 1.602176565e-19 / 10.**30.
    b1w = data_w['BP']

    v0f = data_f['V0']
    b0f = data_f['B0'] * 10.**9. / 1.602176565e-19 / 10.**30.
    b1f = data_f['BP']

    vref = 30.
    bref = 100. * 10.**9. / 1.602176565e-19 / 10.**30.

    if useasymm:
        Vi = 0.94 * v0w
        Vf = 1.06 * v0w
    else:
        Vi = 0.94 * (v0w + v0f) / 2.
        Vf = 1.06 * (v0w + v0f) / 2.

    a3f = 9. * v0f**3. * b0f / 16. * (b1f - 4.)
    a2f = 9. * v0f**(7./3.) * b0f / 16. * (14. - 3. * b1f)
    a1f = 9. * v0f**(5./3.) * b0f / 16. * (3. * b1f - 16.)
    a0f = 9. * v0f * b0f / 16. * (6. - b1f)

    a3w = 9. * v0w**3. * b0w / 16. * (b1w - 4.)
    a2w = 9. * v0w**(7./3.) * b0w / 16. * (14. - 3. * b1w)
    a1w = 9. * v0w**(5./3.) * b0w / 16. * (3. * b1w - 16.)
    a0w = 9. * v0w * b0w / 16. * (6. - b1w)

    x = [0, 0, 0, 0, 0, 0, 0]

    x[0] = (a0f - a0w)**2
    x[1] = 6. * (a1f - a1w) * (a0f - a0w)
    x[2] = -3. * (2. * (a2f - a2w) * (a0f - a0w) + (a1f - a1w)**2.)
    x[3] = -2. * (a3f - a3w) * (a0f - a0w) - 2. * (a2f - a2w) * (a1f - a1w)
    x[4] = -3./5. * (2. * (a3f - a3w) * (a1f - a1w) + (a2f - a2w)**2.)
    x[5] = -6./7. * (a3f - a3w) * (a2f - a2w)
    x[6] = -1./3. * (a3f - a3w)**2.

    y = [0, 0, 0, 0, 0, 0, 0]

    y[0] = (a0f + a0w)**2 / 4.
    y[1] = 3. * (a1f + a1w) * (a0f + a0w) / 2.
    y[2] = -3. * (2. * (a2f + a2w) * (a0f + a0w) + (a1f + a1w)**2.) / 4.
    y[3] = -(a3f + a3w) * (a0f + a0w) / 2. - (a2f + a2w) * (a1f + a1w) / 2.
    y[4] = -3./20. * (2. * (a3f + a3w) * (a1f + a1w) + (a2f + a2w)**2.)
    y[5] = -3./14. * (a3f + a3w) * (a2f + a2w)
    y[6] = -1./12. * (a3f + a3w)**2.

    Fi = np.zeros_like(Vi)
    Ff = np.zeros_like(Vf)

    Gi = np.zeros_like(Vi)
    Gf = np.zeros_like(Vf)

    for n in range(7):
        Fi = Fi + x[n] * Vi**(-(2.*n-3.)/3.)
        Ff = Ff + x[n] * Vf**(-(2.*n-3.)/3.)

        Gi = Gi + y[n] * Vi**(-(2.*n-3.)/3.)
        Gf = Gf + y[n] * Vf**(-(2.*n-3.)/3.)

    Delta = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi))
    Deltarel = 100. * np.sqrt((Ff - Fi) / (Gf - Gi))
    if useasymm:
        Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
                 / v0w / b0w * vref * bref
    else: 
        Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
                 / (v0w + v0f) / (b0w + b0f) * 4. * vref * bref

    return Delta[0], Deltarel[0], Delta1[0]

def get_volumes(n_vol, calc, alat=None):
    abc = np.array(calc["vectors"])
    vol = np.linalg.det(abc) * (calc["alat"] ** 3)
    k = (n_vol - 1) / 2
    return np.linspace(1-0.02*k, 1+0.02*k, n_vol) * vol

def get_alats(volumes, calc):
    vol_abc = np.linalg.det(np.array(calc["vectors"]))
    return (volumes / vol_abc) ** (1./3)


class DeltaCalculation(object):
    
    def __init__(self, settings, uuid, logger=None):
        """ A class for delta factor calculation
        """
        if logger is not None:
            self._log = True
            self._logger = logger
        self._cwd = os.getcwd()
        self.settings = settings
        self.element = settings.calc["element"]
        self.pseudo_file = None
        if self._log:
            self._logger.info("Uuid: {}".format(uuid))
        self._calc_dir = os.path.join(self.element, uuid)
        if not os.path.exists(self._calc_dir):
            os.makedirs(self._calc_dir)
        self._switch_dir()


    def _switch_dir(self):
        if os.getcwd() == self._cwd:
            os.chdir(self._calc_dir)
        elif os.getcwd() == self._calc_dir:
            os.chdir(self._cwd)
        else:
            print "DeltaCalculation._switch_dir: found myself in a strange directory: {}".format(os.getcwd())

    def add_pseudo(self, pseudo_file):
        self.pseudo_file = pseudo_file

    def run_calcs(self, fdf_file):
        volumes = get_volumes(self.settings.volumes, self.settings.calc)
        alats = get_alats(volumes, self.settings.calc)
        x, y = [], []
        siesta_calc = SiestaCalculation(self.settings, self.pseudo_file, fdf_file=fdf_file)
        for alat in alats:
            siesta_calc.prepare(alat)
            if not siesta_calc.is_run:
                siesta_calc.run()
            e = siesta_calc.results()
            if e is not None:
                x.append(float(e[0]))
                y.append(e[1])
        self.volumes = np.array(x) / self.settings.calc["nat"]
        self.energies = np.array(y) / self.settings.calc["nat"]
        if self._log:
            self._logger.debug("Volumes per atom: {}".format(self.volumes))
            self._logger.debug("Energies per atom: {}".format(self.energies))

    def get_delta(self):
        if hasattr(self.settings, 'reference_file'):
            ref_file_name = self.settings.reference_file
        else:
            ref_file_name = "delta/WIEN2k.txt"

        ref_data = read_ref_data(ref_file_name)
        ref_data_el = ref_data[ref_data['element'] == self.element]

        if len(self.volumes) == 7:
            x_p = self.volumes
            y_p = self.energies
        else:
            p = np.polyfit(self.volumes, self.energies, 2)
            min_p = -p[1] / (2*p[0])
            x_p = np.linspace(0.94*min_p, 1.06*min_p, 7)
            y_p = np.polyval(p, x_p)
        vol, bulk_mod, bulk_deriv, _ = BM(np.vstack((x_p, y_p)).T)
        if self._log:
            self._logger.debug("Equil. vol\tBulk modulus\tBulk mod. deriv")
            self._logger.debug("{}\t{}\t{}".format(vol, bulk_mod, bulk_deriv))
        our_data = np.core.records.fromrecords([(self.element, vol, bulk_mod, bulk_deriv), ], names=('element', 'V0', 'B0', 'BP'))
        delta, rel_delta, _ = calcDelta(our_data, ref_data_el, useasymm=False)
        self.delta = delta
        self.rel_delta = rel_delta 
        if self._log:
            self._logger.info("Equilibrium volume per atom =    {min_p:6.6} A^3".format(min_p=vol))
            self._logger.info("""
                                delta, meV/atom  rel_delta, % 
Delta factor                =    {delta:6.4}       {rel_delta:6.4}""".format(delta=delta, rel_delta=rel_delta))
    