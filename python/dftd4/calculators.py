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
"""ASE Calculator implementation for the dftd4 program."""

from __future__ import print_function
from typing import List

from ctypes import c_int, c_double, c_bool

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree, Bohr

import numpy as np

from dftd4.interface import DFTD4Library, P_MBD_APPROX_ATM, P_REFQ_EEQ
__all__ = ["DispersionCorrection", "D4_model"]


class DispersionCorrection(Calculator):
    """Base calculator wrapping a shared library in ctypes."""
    implemented_properties = []  # type: List[str]
    default_parameters = {}  # type: dict

    calc = None
    _debug = False

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False, label=None,
                 atoms=None, library=None, xc=None, calc=None, **kwargs):
        """Construct the base calculator object."""

        Calculator.__init__(self, restart, ignore_bad_restart_file, label, atoms,
                            **kwargs)

        if library is not None:
            self.library = library
        else:
            self.library = DFTD4Library()

        def missing_library_call(**kwargs) -> None:
            """Raise an error if not set by child class."""
            raise NotImplementedError("Child class must replace function call!")

        self.library_call = missing_library_call

        # loads the default parameters and updates with actual values
        self.parameters = self.get_default_parameters()
        # load damping parameters for functional first, than override with set
        if xc is not None:
            self.load_damping_parameters(xc)
        # now set all parameters
        self.set(**kwargs)

        # used as dispersion correction
        if calc is not None:
            self.calc = calc

    def output_file_name(self) -> str:
        """create output file name from label as ctype"""
        if self.label is not None:
            return self.label + ".out"
        return "-"  # standard output

    def get_property(self, name, atoms=None, allow_calculation=True):
        """Taken from the dftd3-calculator in ASE.

        Here we are performing the dispersion correction first,
        if it fails the time for the DFT calculation is not wasted."""

        self_result = Calculator.get_property(self, name, atoms, allow_calculation)
        calc_result = None
        if self.calc is not None:
            calc_result = self.calc.get_property(name, atoms, allow_calculation)
            # little hack for get_potential_energy bypassing this routine
            # in case of force_consistent=True
            if 'free_energy' in self.calc.results:
                # in case we request 'free_energy' with get_property,
                # it will work the first time, the second call will pile
                # up 'free_energy' in the results-array of the dispersion
                # correction and return a wrong value...
                self.results['free_energy'] = self.results['energy'] \
                                              + self.calc.results['free_energy']
                # so here is a really ugly fix to also account for this case
                # -> fix the calculator interface, pls.
                if name == 'free_energy':
                    self_result = self.results['free_energy']
                    calc_result = None
            else:
                del self.results['free_energy']

        if calc_result is None and self_result is None:
            return None
        if calc_result is None:
            return self_result
        if self_result is None:
            return calc_result
        return calc_result + self_result

    def load_damping_parameters(self, functional: str) -> None:
        """create a list of arguments."""
        raise NotImplementedError("Child class must implement damping parameters!")

    def create_arguments(self) -> dict:
        """create a list of arguments."""
        raise NotImplementedError("Child class must implement argument list!")

    def store_results(self, results: dict) -> None:
        """store results after successful calculation."""
        raise NotImplementedError("Child class must implement result processing!")

    # pylint: disable=dangerous-default-value
    def calculate(self, atoms=None, properties: List[str] = None,
                  system_changes: List[str] = all_changes) -> None:
        """calculation interface to libdftd4"""

        if not properties:
            properties = ['energy']
        Calculator.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        kwargs = self.create_arguments()
        results = self.library_call(**kwargs)
        self.store_results(results)


class D4_model(DispersionCorrection):  # pylint: disable=invalid-name
    """actual implementation of the D4 model."""
    implemented_properties = [
        'energy',
        'free_energy',
        'forces',
        'stress',
    ]
    default_parameters = {
        'parallel': 0,
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

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False, label='dftd4',
                 atoms=None, library=None, xc=None, calc=None, **kwargs):
        """Construct the dftd-calculator object."""

        DispersionCorrection.__init__(self, restart, ignore_bad_restart_file,
                                      label, atoms, library, xc, calc, **kwargs)

        self.library_call = self.library.D4Calculation

    def load_damping_parameters(self, functional: str) -> None:
        """load damping parameters and connect it to class"""
        params = self.library.D4Parameters(functional)
        if not params:
            raise RuntimeError("D4 is not parametrized for '{}'"
                               .format(functional))
        self.set(**params)

    def create_arguments(self) -> dict:
        """create a list of arguments."""
        kwargs = {
            'natoms': len(self.atoms),
            'numbers': np.array(self.atoms.get_atomic_numbers(), dtype=c_int),
            'charge': self.atoms.get_initial_charges().sum().round(),
            'positions': np.array(self.atoms.get_positions()/Bohr, dtype=c_double),
            'cell': np.array(self.atoms.get_cell()/Bohr, dtype=c_double),
            'pbc': np.array(self.atoms.get_pbc(), dtype=c_bool),
            'options': self.parameters,
            'output': self.output_file_name(),
        }
        return kwargs

    def store_results(self, results: dict) -> None:
        """store results after successful calculation."""
        self.results['energy'] = results['energy']*Hartree
        self.results['free_energy'] = self.results['energy']
        self.results['forces'] = -results['gradient']*Hartree/Bohr
        if 'stress tensor' in results:
            self.results['stress'] = results['stress tensor']*Hartree/Bohr**3
        self.results['polarizibilities'] = results['polarizibilities']*Bohr**3
        self.results['c6_coefficients'] = results['c6_coefficients']*Hartree*Bohr**6
        self.results['charges'] = results['charges']


class D3_model(DispersionCorrection):  # pylint: disable=invalid-name
    """actual implementation of a D3 like model."""
    implemented_properties = [
        'energy',
        'free_energy',
        'forces',
        'stress',
    ]
    default_parameters = {
        'parallel': 0,
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
        # switches for the shared library, the respective property is simply
        # not calculated if set to False (despite what implemented_properties says)
        'properties': True,
        'energy': True,
        'forces': True,
        'hessian': False,  # unused
    }

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False, label='dftd4',
                 atoms=None, library=None, xc=None, calc=None, **kwargs):
        """Construct the dftd-calculator object."""

        DispersionCorrection.__init__(self, restart, ignore_bad_restart_file,
                                      label, atoms, library, xc, calc, **kwargs)

        self.library_call = self.library.D3Calculation

    def load_damping_parameters(self, functional: str) -> None:
        """Note that there are no official damping parameters for this
        flavour of D3 model, so we will let it pass silently."""

    def create_arguments(self) -> dict:
        """create a list of arguments."""
        kwargs = {
            'natoms': len(self.atoms),
            'numbers': np.array(self.atoms.get_atomic_numbers(), dtype=c_int),
            'positions': np.array(self.atoms.get_positions()/Bohr, dtype=c_double),
            'cell': np.array(self.atoms.get_cell()/Bohr, dtype=c_double),
            'pbc': np.array(self.atoms.get_pbc(), dtype=c_bool),
            'options': self.parameters,
            'output': self.output_file_name(),
        }
        return kwargs

    def store_results(self, results: dict) -> None:
        """store results after successful calculation."""
        self.results['energy'] = results['energy']*Hartree
        self.results['free_energy'] = self.results['energy']
        self.results['forces'] = -results['gradient']*Hartree/Bohr
        if 'stress tensor' in results:
            self.results['stress'] = results['stress tensor']*Hartree/Bohr**3
        self.results['polarizibilities'] = results['polarizibilities']*Bohr**3
        self.results['c6_coefficients'] = results['c6_coefficients']*Hartree*Bohr**6
