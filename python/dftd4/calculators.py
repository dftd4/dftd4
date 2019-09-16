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
nullptr = c_double_p()
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

    def to_dict(self) -> dict:
        """return structure as dictionary"""
        return {'s6': self.s6,
                's8': self.s8,
                's10': self.s10,
                'a1': self.a1,
                'a2': self.a2,
                's9': self.s9,
                'alp': self.alp,
                'beta': self.beta}

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

def load_dftd4_library(library_path: Optional[str] = None) -> CDLL:
    """loader for dftd4 shared object"""
    if library_path is not None:
        lib = cdll.LoadLibrary(library_path)
    else:
        lib = cdll.LoadLibrary('libdftd4.so')

    lib.D4_calculation.argtypes = [
        c_int_p,  # number of atoms
        c_int_p,  # atomic numbers, dimension(number of atoms)
        c_double_p,  # molecular charge
        c_double_p,  # cartesian coordinates, dimension(3*number of atoms)
        c_char_p,  # output file name
        POINTER(DFTD_parameters),
        POINTER(DFTD_options),
        c_double_p,  # energy
        c_double_p,  # gradient, dimension(3*number of atoms)
        c_double_p,  # hessian, dimension((3*number of atoms)**2)
        c_double_p,  # polarizibilities, dimension(number of atoms)
        c_double_p,  # C6-coefficients, dimension(number of atoms**2)
        c_double_p]  # partial charges, dimension(number of atoms)

    lib.D4_PBC_calculation.argtypes = [
        c_int_p,  # number of atoms
        c_int_p,  # atomic numbers, dimension(number of atoms)
        c_double_p,  # molecular charge
        c_double_p,  # cartesian coordinates, dimension(3*number of atoms)
        c_double_p,  # lattice parameters, dimension(9)
        c_bool_p,  # periodicity of the system
        c_char_p,  # output file name
        POINTER(DFTD_parameters),
        POINTER(DFTD_options),
        c_double_p,  # energy
        c_double_p,  # gradient, dimension(3*number of atoms)
        c_double_p,  # lattice gradient, dimension(9)
        c_double_p,  # stress tensor, dimension(6)
        c_double_p,  # hessian, dimension((3*number of atoms)**2)
        c_double_p,  # polarizibilities, dimension(number of atoms)
        c_double_p,  # C6-coefficients, dimension(number of atoms**2)
        c_double_p]  # partial charges, dimension(number of atoms)

    lib.D4_damping_parameters.argtypes = [
        c_char_p,
        POINTER(DFTD_parameters),
        c_int_p]

    lib.D3_calculation.argtypes = [
        c_int_p,  # number of atoms
        c_int_p,  # atomic numbers, dimension(number of atoms)
        c_double_p,  # cartesian coordinates, dimension(3*number of atoms)
        c_char_p,  # output file name
        POINTER(DFTD_parameters),
        POINTER(DFTD_options),
        c_double_p,  # energy
        c_double_p,  # gradient, dimension(3*number of atoms)
        c_double_p,  # hessian, dimension((3*number of atoms)**2)
        c_double_p,  # polarizibilities, dimension(number of atoms)
        c_double_p]  # C6-coefficients, dimension(number of atoms**2)

    lib.D3_PBC_calculation.argtypes = [
        c_int_p,  # number of atoms
        c_int_p,  # atomic numbers, dimension(number of atoms)
        c_double_p,  # cartesian coordinates, dimension(3*number of atoms)
        c_double_p,  # lattice parameters, dimension(9)
        c_bool_p,  # periodicity of the system
        c_char_p,  # output file name
        POINTER(DFTD_parameters),
        POINTER(DFTD_options),
        c_double_p,  # energy
        c_double_p,  # gradient, dimension(3*number of atoms)
        c_double_p,  # lattice gradient, dimension(9)
        c_double_p,  # stress tensor, dimension(6)
        c_double_p,  # hessian, dimension((3*number of atoms)**2)
        c_double_p,  # polarizibilities, dimension(number of atoms)
        c_double_p]  # C6-coefficients, dimension(number of atoms**2)

    return lib

def load_damping_parameters(lib: CDLL, functional: str,
                            lmbd: int = P_MBD_APPROX_ATM) -> dict:
    """load damping parameters from dftd4 shared library"""
    dparam = DFTD_parameters()
    stat: c_int = lib.D4_damping_parameters(functional.encode('utf-8'),
                                            dparam, c_int(lmbd))
    if stat == 0:
        return dparam.to_dict()
    return {}

class D4_model(Calculator):  # pylint: disable=invalid-name
    """Base-class for dftd4 calculators"""
    implemented_properties: List[str] = [
        'energy', 'free_energy', 'forces', 'stress']
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
        'forces': True,
        'charges': False,
        'polarizibilities': False,
        'c6_coefficients': False,
    }

    _debug: bool = False

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=None, atoms=None, **kwargs):
        """Construct the dftd-calculator object."""

        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms, **kwargs)

        # loads the default parameters and updates with actual values
        self.parameters: dict = self.get_default_parameters()
        self.set(**kwargs)

        self._set_implemented_property('forces')
        self._set_implemented_property('charges')
        self._set_implemented_property('polarizibilities')
        self._set_implemented_property('c6_coefficients')

        self._find_lib()
        self._lib: CDLL = load_dftd4_library(self.parameters['libdftd4_path'])

    def _set_implemented_property(self, name) -> None:
        """add properties to implemented properties list"""
        if name in self.parameters:
            if name not in self.implemented_properties and self.parameters[name]:
                self.implemented_properties.append(name)
            if name in self.implemented_properties and not self.parameters[name]:
                self.implemented_properties.remove(name)

    def _find_lib(self) -> None:
        """find libdftd4 if not further specified"""
        # load dftd-library
        if self.parameters['libdftd4_path'] is None:
            from ctypes.util import find_library
            self.parameters['libdftd4_path'] = find_library('dftd4')

    def load_damping_parameters(self, functional) -> None:
        """load damping parameters and connect it to class"""
        params = load_damping_parameters(self._lib, functional)
        if not params:
            raise RuntimeError("D4 is not parametrized for '"+functional+"'")
        self.set(**params)

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
                            'energy' in self.implemented_properties,
                            'forces' in self.implemented_properties,
                            False,  # Hessian
                            self.parameters['print_level'])

    def ctypes_dftd_parameters(self) -> DFTD_parameters:
        """create DFTD_parameters from parameters dictionary as ctype"""
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

    def ctypes_array(self, name, size, ctype=c_double_p):
        """return an optional array as ctype"""
        if name in self.implemented_properties:
            array = np.zeros(size)
            array_p = array.ctypes.data_as(ctype)
        else:
            array = None
            array_p = nullptr
        return (array, array_p)

    def store_result(self, name: str, result, convert=1.0) -> None:
        """store results in dict, make sure the property is implemented"""
        if name in self.implemented_properties:
            self.results[name] = result * convert

    # pylint: disable=too-many-locals, dangerous-default-value
    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """calculation interface to libdftd4"""

        if not properties:
            properties = ['energy']
        Calculator.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

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
        gradient, grad_p = self.ctypes_array('forces', (natoms.value, 3))
        gradlatt = np.zeros((3, 3), order='F')
        glat_p = gradlatt.ctypes.data_as(c_double_p)

        charges, charges_p = self.ctypes_array('charges', natoms.value)
        alphas, alphas_p = self.ctypes_array('polarizibilities', natoms.value)
        c6_coeff, c6_coeff_p = self.ctypes_array('c6_coefficients',
                                                 (natoms.value, natoms.value))

        if self.atoms.pbc.any():
            stat = self._lib.D4_PBC_calculation(
                natoms, attyp_p, charge, coord_p, dlat_p, pbc_p, outfile,
                dparam, opt, energy, grad_p, glat_p, nullptr, nullptr,
                alphas_p, c6_coeff_p, charges_p)
        else:
            stat = self._lib.D4_calculation(
                natoms, attyp_p, charge, coord_p, outfile,
                dparam, opt, energy, grad_p, nullptr,
                alphas_p, c6_coeff_p, charges_p)
        # in case Fortran behaves we find a useful return value
        if stat != 0:
            raise RuntimeError("dftd4 terminated in error.")

        self.store_result('energy', energy.value, Hartree)
        self.store_result('free_energy', energy.value, Hartree)
        self.store_result('forces', gradient, - Hartree / Bohr)
        if self.atoms.pbc.any():
            # dftd4 returns the cell gradient, so we have to backtransform
            # taking into account that the stress tensor should be symmetric,
            # we interface F and C ordered arrays and have to transpose we
            # can simply drop all reshape and transpose steps by writing
            stress = np.dot(self.atoms.cell,
                            gradlatt * Hartree / Bohr / self.atoms.get_volume())
            self.store_result('stress', stress.flat[[0, 4, 8, 5, 2, 1]])

        self.store_result('charges', charges)
        self.store_result('polarizibilities', alphas, Bohr**3)
        self.store_result('c6_coefficients', c6_coeff, Bohr**6 * Hartree)

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
        'forces': True,
        'polarizibilities': False,
        'c6_coefficients': False,
    }

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=None, atoms=None, **kwargs):
        """Construct the dftd-calculator object."""

        D4_model.__init__(self, restart, ignore_bad_restart_file,
                          label, atoms, **kwargs)

    def ctypes_dftd_options(self) -> DFTD_options:
        """create DFTD_options from parameters dictionary as ctype"""
        return DFTD_options(P_MBD_APPROX_ATM,
                            P_REFQ_EEQ,
                            self.parameters['wf'],
                            0.0,   # keep zero for D3-like model
                            0.0,   # keep zero for D3-like model
                            True,  # properties
                            'energy' in self.implemented_properties,
                            'forces' in self.implemented_properties,
                            False,  # Hessian
                            self.parameters['print_level'])

    # pylint: disable=too-many-locals, dangerous-default-value
    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """calculation interface to libdftd4"""

        if not properties:
            properties = ['energy']
        Calculator.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        natoms: c_int = self.ctypes_number_of_atoms()
        attyp_p: pointer = self.ctypes_atomic_numbers()
        coord_p: pointer = self.ctypes_positions()
        dlat_p: pointer = self.ctypes_cell()
        pbc_p: pointer = self.ctypes_pbc()
        opt: DFTD_options = self.ctypes_dftd_options()
        dparam: DFTD_parameters = self.ctypes_dftd_parameters()
        outfile: c_char_p = self.ctypes_output_file()

        energy = c_double(0.0)
        gradient, grad_p = self.ctypes_array('forces', (natoms.value, 3))
        gradlatt = np.zeros((3, 3), order='F')
        glat_p = gradlatt.ctypes.data_as(c_double_p)

        alphas, alphas_p = self.ctypes_array('polarizibilities', natoms.value)
        c6_coeff, c6_coeff_p = self.ctypes_array('c6_coefficients',
                                                 (natoms.value, natoms.value))

        if atoms.pbc.any():
            stat = self._lib.D3_PBC_calculation(
                natoms, attyp_p, coord_p, dlat_p, pbc_p, outfile,
                dparam, opt, energy, grad_p, glat_p, nullptr, nullptr,
                alphas_p, c6_coeff_p)
        else:
            stat = self._lib.D3_calculation(
                natoms, attyp_p, coord_p, outfile,
                dparam, opt, energy, grad_p, nullptr,
                alphas_p, c6_coeff_p)
        # in case Fortran behaves we find a useful return value
        if stat != 0:
            raise RuntimeError("dftd4 terminated in error.")

        self.store_result('energy', energy.value, Hartree)
        self.store_result('free_energy', energy.value, Hartree)
        self.store_result('forces', gradient, - Hartree / Bohr)
        if self.atoms.pbc.any():
            # dftd4 returns the cell gradient, so we have to backtransform
            # taking into account that the stress tensor should be symmetric,
            # we interface F and C ordered arrays and have to transpose we
            # can simply drop all reshape and transpose steps by writing
            stress = np.dot(self.atoms.cell,
                            gradlatt * Hartree / Bohr / self.atoms.get_volume())
            self.results['stress'] = stress.flat[[0, 4, 8, 5, 2, 1]]

        self.store_result('polarizibilities', alphas, Bohr**3)
        self.store_result('c6_coefficients', c6_coeff, Bohr**6 * Hartree)
