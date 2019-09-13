# This file is part of dftd4.
#
# Copyright (C) 2019 Sebastian Ehlert
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.
"""
Wrapper for the dftd4-program.
"""

from __future__ import print_function

from typing import List, Optional

from ctypes import cdll, CDLL, Structure, c_double, c_int, c_bool, \
                   POINTER, pointer, c_char_p

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree, Bohr

import numpy as np

# pylint: disable=invalid-name
c_int_p = POINTER(c_int)
c_bool_p = POINTER(c_bool)
c_double_p = POINTER(c_double)
# pylint: enable=invalid-name

P_REFQ_EEQ: int = 5

P_MBD_NONE: int = 0
P_MBD_RPALIKE: int = 1
P_MBD_APPROX_ATM: int = 3

class DFTD_parameters(Structure): # pylint: disable=invalid-name, too-few-public-methods
    """DFTD damping parameters"""
    _fields_ = [
        ('s6', c_double),
        ('s8', c_double),
        ('s10', c_double),
        ('a1', c_double),
        ('a2', c_double),
        ('s9', c_double),
        ('alp', c_int),
        ('beta', c_double),
    ]

class DFTD_options(Structure): # pylint: disable=invalid-name, too-few-public-methods
    """Options for the dftd4 library"""
    _fields_ = [
        ('lmbd', c_int),
        ('refq', c_int),
        ('wf', c_double),
        ('g_a', c_double),
        ('g_c', c_double),
        ('lmolpol', c_bool),
        ('lenergy', c_bool),
        ('lgradient', c_bool),
        ('lhessian', c_bool),
        ('print_level', c_int),
    ]

class D4_model(Calculator):  # pylint: disable=invalid-name
    """Base-class for dftd4 calculators"""
    implemented_properties: List[str] = ['energy', 'free_energy', 'forces']
    default_parameters = {
        'parallel': 0,
        'gradient': True,
        'libdftd4_path': None,
        's6': 1.0,
        's8': 1.0,
        's10': 0.0,
        's9': 1.0,
        'a1': 0.4,
        'a2': 4.5,
        'alp': 14,
        'wf': 6.0,
        'g_a': 3.0,
        'g_c': 2.0,
        'print_level': 2,
    }

    _lib: Optional[CDLL] = None
    _debug: bool = False

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=None, atoms=None, **kwargs):
        """Construct the dftd-calculator object."""

        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms, **kwargs)

        # loads the default parameters and updates with actual values
        self.parameters: dict = self.get_default_parameters()
        self.set(**kwargs)

    def _find_lib(self) -> None:
        # load dftd-library
        if self.parameters['libdftd4_path'] is None:
            from ctypes.util import find_library
            self.parameters['libdftd4_path'] = find_library('dftd4')

    def _open_lib(self):

        if self.parameters['libdftd4_path'] is not None:
            self._lib: CDLL = cdll.LoadLibrary(self.parameters['libdftd4_path'])
        else:
            self._lib: CDLL = cdll.LoadLibrary('libdftd4.so')

        self._lib.D4_calculation.argtypes = [
            c_int_p, c_int_p, c_double_p, c_double_p, c_char_p,
            POINTER(DFTD_parameters), POINTER(DFTD_options),
            c_double_p, c_double_p]

        self._lib.D4_PBC_calculation.argtypes = [
            c_int_p, c_int_p, c_double_p, c_double_p, c_double_p, c_bool_p,
            c_char_p, POINTER(DFTD_parameters), POINTER(DFTD_options),
            c_double_p, c_double_p, c_double_p]

        self._lib.D3_calculation.argtypes = [
            c_int_p, c_int_p, c_double_p, c_char_p,
            POINTER(DFTD_parameters), POINTER(DFTD_options),
            c_double_p, c_double_p]

        self._lib.D3_PBC_calculation.argtypes = [
            c_int_p, c_int_p, c_double_p, c_double_p, c_bool_p,
            c_char_p, POINTER(DFTD_parameters), POINTER(DFTD_options),
            c_double_p, c_double_p, c_double_p]

    def ctypes_output_file(self) -> c_char_p:
        """create output file name from label as ctype"""
        if self.label is None:
            self.label: str = "saw"
        output_file: str = self.label + ".out"
        return c_char_p(output_file.encode('utf-8'))

    def ctypes_dftd_options(self) -> DFTD_options:
        """create DFTD_options from parameters dictionary as ctype"""
        return DFTD_options(P_MBD_APPROX_ATM,
                            P_REFQ_EEQ,
                            self.parameters['wf'],
                            self.parameters['g_a'],
                            self.parameters['g_c'],
                            True,  # properties
                            True,  # energy
                            True,  # gradient
                            False,
                            self.parameters['print_level'])

    def ctypes_dftd_parameters(self) -> DFTD_parameters:
        """create DFTD_options from parameters dictionary as ctype"""
        return DFTD_parameters(self.parameters['s6'],
                               self.parameters['s8'],
                               self.parameters['s10'],
                               self.parameters['a1'],
                               self.parameters['a2'],
                               self.parameters['s9'],
                               self.parameters['alp'],
                               1.0)

    def ctypes_number_of_atoms(self) -> c_int:
        """return number of atoms from current Atoms-object as ctype"""
        return c_int(len(self.atoms))

    def ctypes_atomic_numbers(self) -> pointer:
        """return atomic numbers from current Atoms-object as ctype"""
        attyp = np.array(self.atoms.get_atomic_numbers(), dtype=np.int32)
        return attyp.ctypes.data_as(c_int_p)

    def ctypes_positions(self) -> pointer:
        """return Cartesian coordinates from current Atoms-object as ctype"""
        coord = self.atoms.get_positions()/Bohr # from Angstrom to Bohr
        return coord.ctypes.data_as(c_double_p)

    def ctypes_cell(self) -> pointer:
        """return lattice parameter from current Atoms-object as ctype"""
        lattice = self.atoms.get_cell()/Bohr # from Angstrom to Bohr
        return lattice.ctypes.data_as(c_double_p)

    def ctypes_initial_charge(self) -> c_double:
        """return total charge of the system from current Atoms-object as ctype"""
        return c_double(self.atoms.get_initial_charges().sum().round())

    def ctypes_pbc(self) -> pointer:
        """returns the periodic directions from current Atoms-object as ctype"""
        pbc = self.atoms.get_pbc()
        return pbc.ctypes.data_as(c_bool_p)

    # pylint: disable=too-many-locals, dangerous-default-value
    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """calculation interface to libdftd4"""

        if not properties:
            properties = ['energy']
        Calculator.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        if self._lib is None:
            self._open_lib()

        natoms: c_int = self.ctypes_number_of_atoms()
        attyp_p: pointer = self.ctypes_atomic_numbers()
        coord_p: pointer = self.ctypes_positions()
        dlat_p: pointer = self.ctypes_cell()
        charge: c_double = self.ctypes_initial_charge()
        pbc_p: pointer = self.ctypes_pbc()
        opt: DFTD_options = self.ctypes_dftd_options()
        dparam: DFTD_parameters = self.ctypes_dftd_parameters()
        outfile: c_char_p = self.ctypes_output_file()

        energy = c_double(0.0)
        gradient = np.zeros((3, natoms.value), order='F')
        grad_p = gradient.ctypes.data_as(c_double_p)
        gradlatt = np.zeros((3, 3), order='F')
        glat_p = gradlatt.ctypes.data_as(c_double_p)

        if atoms.pbc.any():
            stat = self._lib.D4_PBC_calculation(
                natoms, attyp_p, charge, coord_p, dlat_p, pbc_p, outfile,
                dparam, opt, energy, grad_p, glat_p)
        else:
            stat = self._lib.D4_calculation(
                natoms, attyp_p, charge, coord_p, outfile,
                dparam, opt, energy, grad_p)
        # in case Fortran behaves we find a useful return value
        if stat != 0:
            raise RuntimeError("dftd4 terminated in error.")

        self.results['energy'] = energy.value * Hartree
        self.results['free_energy'] = energy.value * Hartree
        self.results['forces'] = - gradient.T * Hartree / Bohr
        if self.atoms.pbc.any():
            # dftd4 returns the cell gradient, so we have to backtransform
            # taking into account that the stress tensor should be symmetric,
            # we interface F and C ordered arrays and have to transpose we
            # can simply drop all reshape and transpose steps by writing
            stress = np.dot(self.atoms.cell,
                            gradlatt * Hartree / Bohr / self.atoms.get_volume())
            self.results['stress'] = stress.flat[[0, 4, 8, 5, 2, 1]]

class D3_model(D4_model):  # pylint: disable=invalid-name
    """D3-like dispersion calculator based on dftd4."""
    default_parameters = {
        'parallel': 0,
        'gradient': True,
        'libdftd4_path': None,
        's6': 1.0,
        's8': 1.0,
        's10': 0.0,
        's9': 1.0,
        'a1': 0.4,
        'a2': 4.5,
        'alp': 14,
        'wf': 4.0,
        'print_level': 2,
    }

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=None, atoms=None, **kwargs):
        """Construct the dftd-calculator object."""

        D4_model.__init__(self, restart, ignore_bad_restart_file,
                          label, atoms, **kwargs)

        # loads the default parameters and updates with actual values
        self.parameters: dict = self.get_default_parameters()
        self.set(**kwargs)

    def ctypes_dftd_options(self) -> DFTD_options:
        """create DFTD_options from parameters dictionary as ctype"""
        return DFTD_options(P_MBD_APPROX_ATM,
                            P_REFQ_EEQ,
                            self.parameters['wf'],
                            0.0,   # keep zero for D3-like model
                            0.0,   # keep zero for D3-like model
                            True,  # properties
                            True,  # energy
                            True,  # gradient
                            False,
                            self.parameters['print_level'])

    # pylint: disable=too-many-locals, dangerous-default-value
    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """calculation interface to libdftd4"""

        if not properties:
            properties = ['energy']
        Calculator.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        if self._lib is None:
            self._open_lib()

        natoms: c_int = self.ctypes_number_of_atoms()
        attyp_p: pointer = self.ctypes_atomic_numbers()
        coord_p: pointer = self.ctypes_positions()
        dlat_p: pointer = self.ctypes_cell()
        pbc_p: pointer = self.ctypes_pbc()
        opt: DFTD_options = self.ctypes_dftd_options()
        dparam: DFTD_parameters = self.ctypes_dftd_parameters()
        outfile: c_char_p = self.ctypes_output_file()

        energy = c_double(0.0)
        gradient = np.zeros((3, natoms.value), order='F')
        grad_p = gradient.ctypes.data_as(c_double_p)
        gradlatt = np.zeros((3, 3), order='F')
        glat_p = gradlatt.ctypes.data_as(c_double_p)

        if atoms.pbc.any():
            stat = self._lib.D3_PBC_calculation(
                natoms, attyp_p, coord_p, dlat_p, pbc_p, outfile,
                dparam, opt, energy, grad_p, glat_p)
        else:
            stat = self._lib.D3_calculation(
                natoms, attyp_p, coord_p, outfile,
                dparam, opt, energy, grad_p)
        # in case Fortran behaves we find a useful return value
        if stat != 0:
            raise RuntimeError("dftd4 terminated in error.")

        self.results['energy'] = energy.value * Hartree
        self.results['free_energy'] = energy.value * Hartree
        self.results['forces'] = - gradient.T * Hartree / Bohr
        if self.atoms.pbc.any():
            # dftd4 returns the cell gradient, so we have to backtransform
            # taking into account that the stress tensor should be symmetric,
            # we interface F and C ordered arrays and have to transpose we
            # can simply drop all reshape and transpose steps by writing
            stress = np.dot(self.atoms.cell,
                            gradlatt * Hartree / Bohr / self.atoms.get_volume())
            self.results['stress'] = stress.flat[[0, 4, 8, 5, 2, 1]]
