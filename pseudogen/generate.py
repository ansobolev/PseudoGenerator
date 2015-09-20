#!/usr/bin/env python

"""
generate.py contains classes and functions needed to generate and test pseudopotential
"""


import os
import shutil
import subprocess
import numpy as np
from scipy.optimize import minimize
from settings import *

orbitals = [(1, 0),(2, 0),(2, 1),
            (3, 0),(3, 1),(4, 0),
            (3, 2),(4, 1),(5, 0),
            (4, 2),(5, 1),(6, 0),
            (4, 3),(5, 2),(6, 1),
            (7, 0),(5, 3),(6, 2),
            ]


class InputFile(object):

    _execute_script = None

    def __init__(self, element, n_core, n_val, **kwds):
        self.element = element
        self.xc = kwds.get("xc", "ca")
        self.n_core = n_core
        self.n_val = n_val
        self._is_spin_pol = kwds.get("is_spin_pol", False)
        self._calc_dir = ".".join((self.element, self.xc, self.calc_type))
        self._file_name = ".".join((self._calc_dir, "inp"))
        self._electrons = []

    def set_calc_dir(self, cdir):
        self._calc_dir = os.path.join(cdir, self._calc_dir)
        self._file_name = ".".join((self._calc_dir, "inp"))
        return self._calc_dir

    def add_electrons(self, electrons):
        n_val = len(self._electrons)
        if n_val == self.n_val:
            print self.__class__.__name__ + ".add_electrons: no valence orbitals to be filled!"
            return False
        self._electrons.append(electrons)
        return True

    def _add_electron_lines(self):
        # adding electrons
        n_val = len(self._electrons)
        if n_val != self.n_val:
            print self.__class__.__name__ + ".execute: %i more valence orbital " \
                                            "configurations must be added" % (self.n_val-n_val)
            return False

        for i_val, el in enumerate(self._electrons):
            n, l = orbitals[self.n_core + i_val]
            if not self._is_spin_pol:
                electrons_str = "{n:5}{l:5}{el:10.3}\n"
                self._input_str += electrons_str.format(n=n, l=l, el=float(el))
            else:
                electrons_str = "{n:5}{l:5}{el_up:10.3}{el_dn:10.3}\n"
                el_up, el_dn = [float(f) for f in el]
                self._input_str += electrons_str.format(n=n, l=l, el_up=el_up, el_dn=el_dn)
        return True

    def add_radii(self, *args):
        pass

    def _add_radii_lines(self):
        return True

    def add_lines(self):
        if not (self._add_electron_lines() and self._add_radii_lines()):
            return False
        return True

    def _pre_execute(self):
        if not self.add_lines():
            return False

        if os.path.exists(self._calc_dir):
            shutil.rmtree(self._calc_dir)
            print self.__class__.__name__ + ".execute: removed old calculation directory"

        os.environ['ATOM_PROGRAM'] = ATOM_PROGRAM
        os.environ['ATOM_UTILS_DIR'] = ATOM_UTILS_DIR
        with open(self._file_name, "w") as f:
            f.write(str(self))

    def execute(self):
        self._pre_execute()
        subprocess.check_call([self._execute_script, self._file_name], env=os.environ)
        return self._post_execute()

    def _post_execute(self):
        pass


class AEInputFile(InputFile):

    _execute_script = os.path.join(ATOM_UTILS_DIR, "ae.sh")
    
    def __init__(self, element, n_core, n_val, **kwds):
        self._input_str = ("#\n"
                           "   {calc_type:<2} {title:<49}\n"
                           "   {elem:<2}   {xc:<2}{xc_opt:<1}\n"
                           "         0\n"
                           "{n_core:5}{n_val:5}\n"
                          )
        self.calc_type = "ae"
        self.title = kwds.get("title", element + " all-electron")
        super(AEInputFile, self).__init__(element, n_core, n_val, **kwds)

    def __str__(self):
        return self._input_str.format(calc_type = self.calc_type,
                                      title = self.title,
                                      elem = self.element,
                                      xc = self.xc,
                                      xc_opt = "r" if self._is_spin_pol else " ",
                                      n_core = self.n_core,
                                      n_val = self.n_val)

class PGInputFile(InputFile):

    _execute_script = os.path.join(ATOM_UTILS_DIR, "pg.sh")

    def __init__(self, element, n_core, n_val, core=False, **kwds):
        self._input_str = ("#\n"
                           "   {calc_type:<2} {title:<49}\n"
                           "        {ps_flavor:<3}{ps_radius:9.3}\n"
                           "   {elem:<2}   {xc:<2}{xc_opt:<1}\n"
                           "         0\n"
                           "{n_core:5}    4\n"
                          )
        
        self.calc_type = "pe" if core else "pg"
        self.title = kwds.get("title", element + " pseudopotential generation")
        self.ps_flavor = kwds.get("ps_flavor", "tm2")
        self.ps_radius = kwds.get("ps_radius", 2.0)
        self._radii = None
        super(PGInputFile, self).__init__(element, n_core, n_val, **kwds)

    @property
    def siesta_pp_file(self):
        return ".".join((self.element, self.xc, self.calc_type, "psf"))

    def _add_electron_lines(self):
        # adding electrons
        n_val = len(self._electrons)
        if n_val != self.n_val:
            print self.__class__.__name__ + ".execute: %i more valence orbital " \
                                            "configurations must be added" % (self.n_val-n_val)
            return False
        # adding empty pseudopotential orbitals
        self._electrons += [(0, 0) if self._is_spin_pol else 0 for _ in range(4-n_val)]
        ps_orbitals = orbitals[self.n_core:self.n_core+n_val]
        ls = [orb[1] for orb in ps_orbitals]
        i = 0
        while len(ps_orbitals) < 4:
            orb = orbitals[self.n_core+n_val+i]
            if orb[1] not in ls:
                ls.append(orb[1])
                ps_orbitals.append(orb)
            i += 1

        for i_val, el in enumerate(self._electrons):
            n, l = ps_orbitals[i_val]
            if not self._is_spin_pol:
                electrons_str = "{n:5}{l:5}{el:10.3}\n"
                self._input_str += electrons_str.format(n=n, l=l, el=float(el))
            else:
                electrons_str = "{n:5}{l:5}{el_up:10.3}{el_dn:10.3}\n"
                el_up, el_dn = [float(f) for f in el]
                self._input_str += electrons_str.format(n=n, l=l, el_up=el_up, el_dn=el_dn)
        return True


    def add_radii(self, r_s, r_p, r_d, r_f, r_ps):
        self._radii = [r_s, r_p, r_d, r_f, r_ps]
    
    def _add_radii_lines(self):
        if self._radii is None:
            print self.__class__.__name__ + ".execute: pseudopotential radii must be added"
            return False
        r_s, r_p, r_d, r_f, r_ps = self._radii 

        radii_str = "{r_s:10.5}{r_p:10.5}{r_d:10.5}{r_f:10.5}{flag:10.5}{r_ps:10.5}\n"
        self._input_str += radii_str.format(r_s = r_s,
                                            r_p = r_p,
                                            r_d = r_d,
                                            r_f = r_f,
                                            flag = 0.,
                                            r_ps = r_ps)
        return True

    def _post_execute(self):
        with open(os.path.join(self._calc_dir, 'OUT'), 'r') as f:
            out_lines = [l for l in f.readlines() if '&v' in l]
        ae_flag = True
        ae = []
        ps = []
        data_type = [("nl", "|S2"),
                     ("s", np.float),
                     ("occ", np.float),
                     ("eigenvalue", np.float),
                     ("kin", np.float),
                     ("pot", np.float)]
        for line in out_lines[1:]:
            if '----' in line:
                ae_flag = False
                continue
            if ae_flag:
                ae.append(tuple(line.split()[:-1]))
            else:
                ps.append(tuple(line.split()[:-1]))
        self.ae = np.core.records.array(ae, dtype=data_type)
        self.ps = np.core.records.array(ps, dtype=data_type)
        err = self.ae['eigenvalue'] - self.ps['eigenvalue']
        return np.sum(err * err) ** 0.5

    def __str__(self):
        return self._input_str.format(calc_type = self.calc_type,
                                      title = self.title,
                                      ps_flavor = self.ps_flavor,
                                      ps_radius = self.ps_radius,
                                      elem = self.element,
                                      xc = self.xc,
                                      xc_opt = "r" if self._is_spin_pol else " ",
                                      n_core = self.n_core,
                                      n_val = self.n_val)

class PTInputFile(InputFile):

    _execute_script = os.path.join(ATOM_UTILS_DIR, "pt.sh")

    def __init__(self, element, n_core, n_val, **kwds):
        self.title = kwds.get("title", element + " pseudopotential test")
        self.calc_type = "pt"
        super(PTInputFile, self).__init__(element, n_core, n_val, **kwds)
        pp_calc_type = "pe" if kwds.get("core", False) else "pg"
        self._pp_file_name = ".".join((self.element, self.xc, pp_calc_type, "vps"))
        assert os.path.exists(self._pp_file_name)
        self._calc_dir = self._calc_dir + "-" + ".".join((self.element, self.xc, pp_calc_type))
        self._configurations = []

    def add_configuration(self, electrons):
        conf = AEInputFile(self.element, 
                           self.n_core,
                           len(electrons),
                           xc=self.xc,
                           is_spin_pol=self._is_spin_pol)
        for el in electrons:
            conf.add_electrons(el)
        self._configurations.append(conf)

    def _pre_execute(self):
        if os.path.exists(self._calc_dir):
            shutil.rmtree(self._calc_dir)
            print self.__class__.__name__ + ".execute: removed old calculation directory"

        os.environ['ATOM_PROGRAM'] = ATOM_PROGRAM
        os.environ['ATOM_UTILS_DIR'] = ATOM_UTILS_DIR
        with open(self._file_name, "w") as f:
            for conf in self._configurations:
                conf.add_lines()
                f.write(str(conf))
                conf.calc_type = self.calc_type
            for conf in self._configurations:
                f.write(str(conf))

    def execute(self):
        self._pre_execute()
        subprocess.check_call([self._execute_script, self._file_name, self._pp_file_name], 
                              env=os.environ)
        return self._post_execute()

    def _post_execute(self):
        with open(os.path.join(self._calc_dir, 'OUT'), 'r') as f:
            out_lines = [l for l in f.readlines() if '&d' in l and '&v' not in l]
        # cross-excitations
        nconfs = len(self._configurations)
        self.xx = np.zeros((2, nconfs, nconfs))
        ps_flag = False
        head_flag = False
        for line in out_lines[2:]:
            if 'total' in line:
                ps_flag = True
                head_flag = True
                continue
            if head_flag:
                head_flag = False
                continue
            data = line.split()[1:]
            idx = int(data[0]) - 1
            for ix, x in enumerate([float(x) for x in data[1:]]):
                self.xx[ps_flag, ix, idx] = x
                self.xx[ps_flag, idx, ix] = x

        err = self.xx[0] - self.xx[1]
        err_mean = np.sum(np.abs(err)) / (nconfs * (nconfs-1))
        err_max = np.max(np.abs(err))
        # return mean value
        # return err_mean
        print 'Mean: %f' % (err_mean)
        print 'Max: %f' % (err_max)
        return err_mean, err_max


def generate_pseudo(calc, electrons, radii):
    ps = PGInputFile(**calc)
    for e in electrons:
        ps.add_electrons(e)
    ps.add_radii(*radii)
    file_name = ps.siesta_pp_file
    return file_name, ps.execute()

def test_pseudo(calc, configs):
    pt = PTInputFile(**calc)
    for c in configs:
        pt.add_configuration(c)
    return pt.execute()    

def generate_test_cycle(*args):
    generate_pseudo(calc, electrons, args)
    return test_pseudo(calc, configs)

