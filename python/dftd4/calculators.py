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
from ctypes.util import find_library

from ase.calculators.calculator import Calculator, all_changes, \
                                       CalculatorSetupError
from ase.units import Hartree, Bohr

import numpy as np

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
        ('properties', c_bool),
        ('energy', c_bool),
        ('forces', c_bool),
        ('hessian', c_bool),
        ('print_level', c_int),
    ]

def load_dftd4_library(library_path: Optional[str] = None) -> CDLL:
    """loader for dftd4 shared object"""
    if library_path is not None:
        lib = cdll.LoadLibrary(library_path)
    else:
        lib = cdll.LoadLibrary('libdftd4.so')

    c_int_p = POINTER(c_int)
    c_bool_p = POINTER(c_bool)
    c_double_p = POINTER(c_double)

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

class SharedObjectCalculator(Calculator):
    """Base class for calculators that wrap an external shared library"""

    library = None

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False, label=None,
                 atoms=None, library=None, **kwargs):
        """Constructor for the shared object based calculator."""

        Calculator.__init__(self, restart, ignore_bad_restart_file, label,
                            atoms, **kwargs)

        if library is not None:
            self.library = library

    def get_number_of_atoms_as(self, ctype=c_int) -> c_int:
        """return number of atoms from current Atoms-object as ctype"""
        return ctype(len(self.atoms))

    def get_atomic_numbers_as(self, ctype=c_int) -> pointer:
        """return atomic numbers from current Atoms-object as ctype"""
        attyp = np.array(self.atoms.get_atomic_numbers(), dtype=ctype)
        return attyp.ctypes.data_as(POINTER(ctype))

    def get_positions_as(self, ctype=c_double, convert=1.0) -> pointer:
        """Cartesian coordinates from current Atoms-object as ctype"""
        coord = np.array(self.atoms.get_positions() * convert, dtype=ctype)
        return coord.ctypes.data_as(POINTER(ctype))

    def get_cell_as(self, ctype=c_double, convert=1.0) -> pointer:
        """lattice parameter from current Atoms-object as ctype"""
        lattice = np.array(self.atoms.get_cell() * convert, dtype=ctype)
        return lattice.ctypes.data_as(POINTER(ctype))

    def get_initial_charge_as(self, ctype=c_double) -> c_double:
        """total charge of the system from current Atoms-object as ctype"""
        return ctype(self.atoms.get_initial_charges().sum().round())

    def get_pbc_as(self, ctype=c_bool) -> pointer:
        """periodic directions from current Atoms-object as ctype"""
        pbc = np.array(self.atoms.get_pbc(), dtype=c_bool)
        return pbc.ctypes.data_as(POINTER(ctype))

    # pylint: disable=protected-access
    def get_struct(self, struct, parameters=None, failsave=False) -> Structure:
        """create a Structure from the parameters dictionary.

        It uses dictionary comprehension to filter the keys in the fields list
        of the Structure class and generate the necessary arguments
        for creating the class on-the-fly."""

        if parameters is None:
            parameters = self.parameters

        if failsave:
            kwargs = {key: parameters[key] for key, _ in struct._fields_
                      if key in parameters}
        else:
            kwargs = {key: parameters[key] for key, _ in struct._fields_}

        return struct(**kwargs)

    # pylint: disable=dangerous-default-value
    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        """Actual calculation."""
        Calculator.calculate(self, atoms, properties, system_changes)

        if self.library is None:
            raise CalculatorSetupError(
                "No shared object bound to {} calculator!"
                .format(self.__class__.__name__))

class D4_model(SharedObjectCalculator):  # pylint: disable=invalid-name
    """Base-class for dftd4 calculators"""
    implemented_properties: List[str] = [
        'energy', 'free_energy', 'forces', 'stress',
        'charges', 'polarizibilities', 'c6_coefficients']
    default_parameters = {
        'parallel': 0,
        'gradient': True,
        'libdftd4_path': None,
        's6': -1.0,
        's8': -1.0,
        's10': 0.0,
        's9': 1.0,
        'a1': -1.0,
        'a2': -1.0,
        'alp': 16,
        'beta': 1.0,  # unused
        'lmbd': P_MBD_APPROX_ATM,
        'refq': P_REFQ_EEQ,
        'wf': 6.0,
        'g_a': 3.0,
        'g_c': 2.0,
        'print_level': 2,
        'properties': True,
        'energy': True,
        'forces': True,
        'hessian': False,
    }

    calc = None
    _debug: bool = False

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=None, atoms=None, library=None, xc=None, calc=None,
                 **kwargs):
        """Construct the dftd-calculator object."""

        SharedObjectCalculator.__init__(self, restart, ignore_bad_restart_file,
                                        label, atoms, library, **kwargs)

        if library is None:
            path = find_library('dftd4')
            self.library: CDLL = load_dftd4_library(path)

        # loads the default parameters and updates with actual values
        self.parameters: dict = self.get_default_parameters()
        self.set(**kwargs)

        if xc is not None:
            self.load_damping_parameters(xc)

        # used as dispersion correction
        if calc is not None:
            self.calc = calc

    def load_damping_parameters(self, functional) -> None:
        """load damping parameters and connect it to class"""
        params = load_damping_parameters(self.library, functional)
        if not params:
            raise RuntimeError("D4 is not parametrized for '"+functional+"'")
        self.set(**params)

    def output_file_name(self) -> c_char_p:
        """create output file name from label as ctype"""
        if self.label is None:
            self.label: str = "saw"
        output_file: str = self.label + ".out"
        return c_char_p(output_file.encode('utf-8'))

    def ctypes_array(self, name, size, ctype=c_double):
        """return an optional array as ctype"""
        if name in self.implemented_properties:
            array = np.zeros(size, dtype=ctype)
            a_ptr = array.ctypes.data_as(POINTER(ctype))
            return (array, a_ptr)
        return (None, None)

    def get_property(self, name, atoms=None, allow_calculation=True):
        """Taken from the dftd3-calculator in ASE.

        Here we are performing the dispersion correction first,
        if it fails the time for the DFT calculation is not wasted."""

        self_result = SharedObjectCalculator.get_property(self, name, atoms,
                                                          allow_calculation)
        calc_result = None
        if self.calc is not None:
            calc_result = self.calc.get_property(name, atoms, allow_calculation)
            # little hack for get_potential_energy bypassing this routine
            # in case of force_consistent=True
            if 'free_energy' in self.calc.results:
                self.results['free_energy'] += self.calc.results['free_energy']
            else:
                del self.results['free_energy']

        if calc_result is None and self_result is None:
            return None
        if calc_result is None:
            return self_result
        if self_result is None:
            return calc_result
        return calc_result + self_result

    def store_result(self, name: str, result, convert=1.0) -> None:
        """store results in dict, make sure the property is implemented"""
        if name in self.implemented_properties:
            self.results[name] = result * convert

    # pylint: disable=too-many-locals, dangerous-default-value
    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """calculation interface to libdftd4"""

        if not properties:
            properties = ['energy']
        SharedObjectCalculator.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        # construct ctypes representations from all intent(in) arguments
        natoms: c_int = self.get_number_of_atoms_as(c_int)
        attyp_p: pointer = self.get_atomic_numbers_as(c_int)
        charge: c_double = self.get_initial_charge_as(c_double)
        coord_p: pointer = self.get_positions_as(c_double, convert=1.0/Bohr)
        dlat_p: pointer = self.get_cell_as(c_double, convert=1.0/Bohr)
        pbc_p: pointer = self.get_pbc_as(c_bool)
        outfile: c_char_p = self.output_file_name()
        dparam: DFTD_parameters = self.get_struct(DFTD_parameters)
        opt: DFTD_options = self.get_struct(DFTD_options)

        # now we need to make space for all intent(out) arguments
        # note: the library will not write to the address if we pass None
        energy = c_double(0.0)
        gradient, grad_p = self.ctypes_array('forces', (natoms.value, 3))
        stress, stress_p = self.ctypes_array('stress', (3, 3))
        charges, charges_p = self.ctypes_array('charges', natoms.value)
        alphas, alphas_p = self.ctypes_array('polarizibilities', natoms.value)
        c6_coeff, c6_coeff_p = self.ctypes_array('c6_coefficients',
                                                 (natoms.value, natoms.value))

        stat = self.library.D4_PBC_calculation(
            natoms, attyp_p, charge, coord_p, dlat_p, pbc_p, outfile,
            dparam, opt, energy, grad_p, None, stress_p, None,
            alphas_p, c6_coeff_p, charges_p)
        # in case Fortran behaves we find a useful return value
        if stat != 0:
            raise RuntimeError("dftd4 terminated in error.")

        self.store_result('energy', energy.value, Hartree)
        self.store_result('free_energy', energy.value, Hartree)
        self.store_result('forces', gradient, - Hartree / Bohr)
        self.store_result('stress', stress, Hartree / Bohr**3)
        self.store_result('charges', charges)
        self.store_result('polarizibilities', alphas, Bohr**3)
        self.store_result('c6_coefficients', c6_coeff, Bohr**6 * Hartree)

class D3_model(D4_model):  # pylint: disable=invalid-name
    """D3-like dispersion calculator based on dftd4."""
    implemented_properties: List[str] = [
        'energy', 'free_energy', 'forces', 'stress',
        'polarizibilities', 'c6_coefficients']
    default_parameters = {
        'parallel': 0,
        'gradient': True,
        'libdftd4_path': None,
        's6': -1.0,
        's8': -1.0,
        's10': 0.0,
        's9': 1.0,
        'a1': -1.0,
        'a2': -1.0,
        'alp': 16,
        'beta': 1.0,  # unused
        'lmbd': P_MBD_APPROX_ATM,
        'refq': P_REFQ_EEQ,
        'wf': 4.0,
        'g_a': 0.0,
        'g_c': 0.0,
        'print_level': 2,
        'properties': False,
        'energy': True,
        'forces': True,
        'hessian': False,
    }

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False, label=None,
                 atoms=None, library=None, xc=None, calc=None, **kwargs):
        """Construct the dftd-calculator object."""

        D4_model.__init__(self, restart, ignore_bad_restart_file,
                          label, atoms, library, xc, calc, **kwargs)

    # pylint: disable=too-many-locals, dangerous-default-value
    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """calculation interface to libdftd4"""

        if not properties:
            properties = ['energy']
        SharedObjectCalculator.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        # construct ctypes representations from all intent(in) arguments
        natoms: c_int = self.get_number_of_atoms_as(c_int)
        attyp_p: pointer = self.get_atomic_numbers_as(c_int)
        coord_p: pointer = self.get_positions_as(c_double, convert=1.0/Bohr)
        dlat_p: pointer = self.get_cell_as(c_double, convert=1.0/Bohr)
        pbc_p: pointer = self.get_pbc_as(c_bool)
        opt: DFTD_options = self.get_struct(DFTD_options)
        dparam: DFTD_parameters = self.get_struct(DFTD_parameters)
        outfile: c_char_p = self.output_file_name()

        # now we need to make space for all intent(out) arguments
        # note: the library will not write to the address if we pass None
        energy = c_double(0.0)
        gradient, grad_p = self.ctypes_array('forces', (natoms.value, 3))
        stress, stress_p = self.ctypes_array('stress', (3, 3))
        alphas, alphas_p = self.ctypes_array('polarizibilities', natoms.value)
        c6_coeff, c6_coeff_p = self.ctypes_array('c6_coefficients',
                                                 (natoms.value, natoms.value))

        stat = self.library.D3_PBC_calculation(
            natoms, attyp_p, coord_p, dlat_p, pbc_p, outfile,
            dparam, opt, energy, grad_p, None, stress_p, None,
            alphas_p, c6_coeff_p)
        # in case Fortran behaves we find a useful return value
        if stat != 0:
            raise RuntimeError("dftd4 terminated in error.")

        self.store_result('energy', energy.value, Hartree)
        self.store_result('free_energy', energy.value, Hartree)
        self.store_result('forces', gradient, - Hartree / Bohr)
        self.store_result('stress', stress, Hartree / Bohr**3)
        self.store_result('polarizibilities', alphas, Bohr**3)
        self.store_result('c6_coefficients', c6_coeff, Bohr**6 * Hartree)
