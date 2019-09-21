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

from typing import List, Optional, Tuple, Any, Union

from ctypes import cdll, CDLL, Structure, c_double, c_int, c_bool, \
                   POINTER, pointer, c_char_p
from ctypes.util import find_library

from ase.calculators.calculator import Calculator, all_changes, \
                                       CalculatorSetupError
from ase.units import Hartree, Bohr

import numpy as np

P_REFQ_EEQ = 5

P_MBD_NONE = 0
P_MBD_RPALIKE = 1
P_MBD_APPROX_ATM = 3


class DFTDParameters(Structure):  # pylint: disable=too-few-public-methods
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


class DFTDOptions(Structure):  # pylint: disable=too-few-public-methods
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

    def to_dict(self) -> dict:
        """return structure as dictionary"""
        return {'lmbd': self.lmbd,
                'refq': self.refq,
                'wf': self.wf,
                'g_a': self.g_a,
                'g_c': self.g_c,
                'properties': self.properties,
                'energy': self.energy,
                'forces': self.forces,
                'hessian': self.hessian,
                'print_level': self.print_level}


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
        POINTER(DFTDParameters),
        POINTER(DFTDOptions),
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
        POINTER(DFTDParameters),
        POINTER(DFTDOptions),
        c_double_p,  # energy
        c_double_p,  # gradient, dimension(3*number of atoms)
        c_double_p,  # lattice gradient, dimension(9)
        c_double_p,  # stress tensor, dimension(9)
        c_double_p,  # hessian, dimension((3*number of atoms)**2)
        c_double_p,  # polarizibilities, dimension(number of atoms)
        c_double_p,  # C6-coefficients, dimension(number of atoms**2)
        c_double_p]  # partial charges, dimension(number of atoms)

    lib.D4_damping_parameters.argtypes = [
        c_char_p,
        POINTER(DFTDParameters),
        c_int_p]

    lib.D3_calculation.argtypes = [
        c_int_p,  # number of atoms
        c_int_p,  # atomic numbers, dimension(number of atoms)
        c_double_p,  # cartesian coordinates, dimension(3*number of atoms)
        c_char_p,  # output file name
        POINTER(DFTDParameters),
        POINTER(DFTDOptions),
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
        POINTER(DFTDParameters),
        POINTER(DFTDOptions),
        c_double_p,  # energy
        c_double_p,  # gradient, dimension(3*number of atoms)
        c_double_p,  # lattice gradient, dimension(9)
        c_double_p,  # stress tensor, dimension(9)
        c_double_p,  # hessian, dimension((3*number of atoms)**2)
        c_double_p,  # polarizibilities, dimension(number of atoms)
        c_double_p]  # C6-coefficients, dimension(number of atoms**2)

    return lib


def load_d4_damping_parameters(lib: CDLL, functional: str,
                               lmbd: int = P_MBD_APPROX_ATM) -> dict:
    """load damping parameters from dftd4 shared library"""
    dparam = DFTDParameters()
    stat = lib.D4_damping_parameters(functional.encode('utf-8'),
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

    # pylint: disable=protected-access
    def get_struct(self, struct, parameters=None, failsave=False):
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
    implemented_properties = [
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
        # switches for the shared library, the respective property is simply
        # not calculated if set to False (despite what implemented_properties says)
        'properties': True,
        'energy': True,
        'forces': True,
        'hessian': False,  # unused
    }
    # check this parameter list to be positive before calling the shared library
    verify = ['s6', 's8', 's9', 's10', 'a1', 'a2', 'alp', 'wf', 'g_a', 'g_c']

    calc = None
    _debug = False

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='dftd4', atoms=None, library=None, xc=None, calc=None,
                 **kwargs):
        """Construct the dftd-calculator object."""

        SharedObjectCalculator.__init__(self, restart, ignore_bad_restart_file,
                                        label, atoms, library, **kwargs)

        if library is None:
            path = find_library('dftd4')
            self.library = load_dftd4_library(path)  # type: CDLL

        # loads the default parameters and updates with actual values
        self.parameters = self.get_default_parameters()
        self.set(**kwargs)

        if xc is not None:
            self.load_damping_parameters(xc)

        # used as dispersion correction
        if calc is not None:
            self.calc = calc

    def load_damping_parameters(self, functional: str) -> None:
        """load damping parameters and connect it to class"""
        params = load_d4_damping_parameters(self.library, functional)
        if not params:
            raise RuntimeError("D4 is not parametrized for '{}'"
                               .format(functional))
        self.set(**params)

    def output_file_name(self) -> c_char_p:
        """create output file name from label as ctype"""
        if self.label is not None:
            output_file = self.label + ".out"
        else:
            output_file = "-"  # standard output
        return c_char_p(output_file.encode('utf-8'))

    def ctypes_array(self, name: str, size: Union[tuple, int], ctype=c_double)\
            -> Tuple[Optional[np.ndarray], Optional[pointer]]:
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
                self.results['free_energy'] = self.results['energy'] \
                                              + self.calc.results['free_energy']
            else:
                del self.results['free_energy']

        if calc_result is None and self_result is None:
            return None
        if calc_result is None:
            return self_result
        if self_result is None:
            return calc_result
        return calc_result + self_result

    def store_result(self, name: str, result: Any, convert=1.0) -> None:
        """store results in dict, make sure the property is implemented"""
        if name in self.implemented_properties:
            self.results[name] = result * convert

    def check_parameters(self):
        """sanity check for the parameters to be verified"""
        check_failed = [key for key in self.verify if self.parameters[key] < 0]
        if check_failed:
            raise CalculatorSetupError('Method parameters {} cannot be negative!'
                                       .format(check_failed))

    # pylint: disable=too-many-locals, dangerous-default-value
    def calculate(self, atoms=None, properties: List[str] = None,
                  system_changes: List[str] = all_changes) -> None:
        """calculation interface to libdftd4"""

        self.check_parameters()

        if not properties:
            properties = ['energy']
        SharedObjectCalculator.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        # construct ctypes representations from all intent(in) arguments
        natoms = c_int(len(self.atoms))
        numbers = np.array(self.atoms.get_atomic_numbers(), dtype=c_int)
        charge = c_double(self.atoms.get_initial_charges().sum().round())
        positions = np.array(self.atoms.get_positions() / Bohr, dtype=c_double)
        cell = np.array(self.atoms.get_cell() / Bohr, dtype=c_double)
        pbc = np.array(self.atoms.get_pbc(), dtype=c_bool)
        outfile = self.output_file_name()
        dparam = self.get_struct(DFTDParameters)
        opt = self.get_struct(DFTDOptions)

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
            natoms,
            numbers.ctypes.data_as(POINTER(c_int)),
            charge,
            positions.ctypes.data_as(POINTER(c_double)),
            cell.ctypes.data_as(POINTER(c_double)),
            pbc.ctypes.data_as(POINTER(c_bool)),
            outfile,
            dparam,
            opt,
            energy,
            grad_p,
            None,
            stress_p,
            None,
            alphas_p,
            c6_coeff_p,
            charges_p)
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
    implemented_properties = [
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
    def __init__(self, restart=None, ignore_bad_restart_file=False, label='dftd3',
                 atoms=None, library=None, xc=None, calc=None, **kwargs):
        """Construct the dftd-calculator object."""

        D4_model.__init__(self, restart, ignore_bad_restart_file,
                          label, atoms, library, xc, calc, **kwargs)

    def load_damping_parameters(self, functional: str) -> None:
        """load damping parameters and connect it to class.

        Note that this D3 method does not have official parameters"""
        raise RuntimeError("This kind of D3 has no published damping parameters")

    # pylint: disable=too-many-locals, dangerous-default-value
    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """calculation interface to libdftd4"""

        self.check_parameters()

        if not properties:
            properties = ['energy']
        SharedObjectCalculator.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        # construct ctypes representations from all intent(in) arguments
        natoms = c_int(len(self.atoms))
        numbers = np.array(self.atoms.get_atomic_numbers(), dtype=c_int)
        positions = np.array(self.atoms.get_positions() / Bohr, dtype=c_double)
        cell = np.array(self.atoms.get_cell() / Bohr, dtype=c_double)
        pbc = np.array(self.atoms.get_pbc(), dtype=c_bool)
        outfile = self.output_file_name()
        dparam = self.get_struct(DFTDParameters)
        opt = self.get_struct(DFTDOptions)

        # now we need to make space for all intent(out) arguments
        # note: the library will not write to the address if we pass None
        energy = c_double(0.0)
        gradient, grad_p = self.ctypes_array('forces', (natoms.value, 3))
        stress, stress_p = self.ctypes_array('stress', (3, 3))
        alphas, alphas_p = self.ctypes_array('polarizibilities', natoms.value)
        c6_coeff, c6_coeff_p = self.ctypes_array('c6_coefficients',
                                                 (natoms.value, natoms.value))

        stat = self.library.D3_PBC_calculation(
            natoms,
            numbers.ctypes.data_as(POINTER(c_int)),
            positions.ctypes.data_as(POINTER(c_double)),
            cell.ctypes.data_as(POINTER(c_double)),
            pbc.ctypes.data_as(POINTER(c_bool)),
            outfile,
            dparam,
            opt,
            energy,
            grad_p,
            None,
            stress_p,
            None,
            alphas_p,
            c6_coeff_p)
        # in case Fortran behaves we find a useful return value
        if stat != 0:
            raise RuntimeError("dftd4 terminated in error.")

        self.store_result('energy', energy.value, Hartree)
        self.store_result('free_energy', energy.value, Hartree)
        self.store_result('forces', gradient, - Hartree / Bohr)
        self.store_result('stress', stress, Hartree / Bohr**3)
        self.store_result('polarizibilities', alphas, Bohr**3)
        self.store_result('c6_coefficients', c6_coeff, Bohr**6 * Hartree)
