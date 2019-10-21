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
"""Wrapper around the C-API of the dftd4 shared library."""

from typing import Optional

from ctypes import Structure, c_int, c_double, c_bool, c_char_p, \
                   POINTER, cdll, CDLL

import os.path as op
import numpy as np
from numpy.distutils.misc_util import get_shared_lib_extension

__all__ = [
    'DFTDOptions',
    'DFTDParameters',
    'DFTD4Library',
    'P_REFQ_EEQ',
    'P_MBD_NONE',
    'P_MBD_RPALIKE',
    'P_MBD_APPROX_ATM',
    'load_library',
]

P_REFQ_EEQ = 5

P_MBD_NONE = 0
P_MBD_RPALIKE = 1
P_MBD_APPROX_ATM = 3

# seems like ctypeslib is not always available
try:
    as_ctype = np.ctypeslib.as_ctypes_type  # pylint:disable=invalid-name
except AttributeError:
    as_ctype = None  # pylint:disable=invalid-name


def load_library(libname: str) -> CDLL:
    """load library cross-platform compatible."""

    if not op.splitext(libname)[1]:
        # Try to load library with platform-specific name
        so_ext = get_shared_lib_extension()
        libname_ext = libname + so_ext
    else:
        libname_ext = libname

    return cdll.LoadLibrary(libname_ext)


class _Structure_(Structure):  # pylint: disable=invalid-name,protected-access
    """patch Structure class to allow returning it arguments as dict."""
    def to_dict(self) -> dict:
        """return structure as dictionary."""
        return {key: getattr(self, key) for key, _ in self._fields_}

    def to_list(self) -> list:
        """return structure as list. Order is the same as in structure."""
        return [getattr(self, key) for key, _ in self._fields_]


class DFTDParameters(_Structure_):
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


class DFTDOptions(_Structure_):
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


def check_ndarray(array: np.ndarray, ctype, size: int, name="array") -> None:
    """check if we got the correct array data"""
    if not isinstance(array, np.ndarray):
        raise ValueError("{} must be of type ndarray".format(name))
    if array.size != size:
        raise ValueError("{} does not have the correct size of {}"
                         .format(name, size))
    if as_ctype is not None:
        if as_ctype(array.dtype) != ctype:
            raise ValueError("{} must be of {} compatible type"
                             .format(name, ctype))


class DFTD4Library:
    """wrapper for the dftd4 shared library"""

    _D4_calculation_ = (
        POINTER(c_int),  # number of atoms
        POINTER(c_int),  # atomic numbers, dimension(number of atoms)
        POINTER(c_double),  # molecular charge
        POINTER(c_double),  # cartesian coordinates, dimension(3*number of atoms)
        c_char_p,  # output file name
        POINTER(DFTDParameters),
        POINTER(DFTDOptions),
        POINTER(c_double),  # energy
        POINTER(c_double),  # gradient, dimension(3*number of atoms)
        POINTER(c_double),  # hessian, dimension((3*number of atoms)**2)
        POINTER(c_double),  # polarizibilities, dimension(number of atoms)
        POINTER(c_double),  # C6-coefficients, dimension(number of atoms**2)
        POINTER(c_double),  # partial charges, dimension(number of atoms)
    )

    _D4_PBC_calculation_ = (
        POINTER(c_int),  # number of atoms
        POINTER(c_int),  # atomic numbers, dimension(number of atoms)
        POINTER(c_double),  # molecular charge
        POINTER(c_double),  # cartesian coordinates, dimension(3*number of atoms)
        POINTER(c_double),  # lattice parameters, dimension(9)
        POINTER(c_bool),  # periodicity of the system
        c_char_p,  # output file name
        POINTER(DFTDParameters),
        POINTER(DFTDOptions),
        POINTER(c_double),  # energy
        POINTER(c_double),  # gradient, dimension(3*number of atoms)
        POINTER(c_double),  # lattice gradient, dimension(9)
        POINTER(c_double),  # stress tensor, dimension(9)
        POINTER(c_double),  # hessian, dimension((3*number of atoms)**2)
        POINTER(c_double),  # polarizibilities, dimension(number of atoms)
        POINTER(c_double),  # C6-coefficients, dimension(number of atoms**2)
        POINTER(c_double),  # partial charges, dimension(number of atoms)
    )

    _D4_damping_parameters_ = (
        c_char_p,
        POINTER(DFTDParameters),
        POINTER(c_int),
    )

    _D3_calculation_ = (
        POINTER(c_int),  # number of atoms
        POINTER(c_int),  # atomic numbers, dimension(number of atoms)
        POINTER(c_double),  # cartesian coordinates, dimension(3*number of atoms)
        c_char_p,  # output file name
        POINTER(DFTDParameters),
        POINTER(DFTDOptions),
        POINTER(c_double),  # energy
        POINTER(c_double),  # gradient, dimension(3*number of atoms)
        POINTER(c_double),  # hessian, dimension((3*number of atoms)**2)
        POINTER(c_double),  # polarizibilities, dimension(number of atoms)
        POINTER(c_double),  # C6-coefficients, dimension(number of atoms**2)
    )

    _D3_PBC_calculation_ = (
        POINTER(c_int),  # number of atoms
        POINTER(c_int),  # atomic numbers, dimension(number of atoms)
        POINTER(c_double),  # cartesian coordinates, dimension(3*number of atoms)
        POINTER(c_double),  # lattice parameters, dimension(9)
        POINTER(c_bool),  # periodicity of the system
        c_char_p,  # output file name
        POINTER(DFTDParameters),
        POINTER(DFTDOptions),
        POINTER(c_double),  # energy
        POINTER(c_double),  # gradient, dimension(3*number of atoms)
        POINTER(c_double),  # lattice gradient, dimension(9)
        POINTER(c_double),  # stress tensor, dimension(9)
        POINTER(c_double),  # hessian, dimension((3*number of atoms)**2)
        POINTER(c_double),  # polarizibilities, dimension(number of atoms)
        POINTER(c_double),  # C6-coefficients, dimension(number of atoms**2)
    )

    def __init__(self, library: Optional[CDLL] = None):
        """construct library from CDLL object."""
        if library is None:
            library = load_library('libdftd4')

        self.library = library
        self._set_argtypes_()

    def _set_argtypes_(self) -> None:
        """define all interfaces."""
        self.library.D4_PBC_calculation.argtypes = self._D4_PBC_calculation_
        self.library.D4_calculation.argtypes = self._D4_calculation_
        self.library.D4_damping_parameters.argtypes = self._D4_damping_parameters_
        self.library.D3_PBC_calculation.argtypes = self._D3_PBC_calculation_
        self.library.D3_calculation.argtypes = self._D3_calculation_

    # pylint: disable=invalid-name, too-many-arguments, too-many-locals
    def D4Calculation(self, natoms: int, numbers, positions, options: dict,
                      charge: float = 0.0, output: str = "-",
                      cell=None, pbc=None) -> dict:
        """wrapper for calling the D4 Calculator from the library."""
        periodic = cell is not None and pbc is not None

        check_ndarray(numbers, c_int, natoms, "numbers")
        check_ndarray(positions, c_double, 3*natoms, "positions")
        if periodic:
            check_ndarray(cell, c_double, 9, "cell")
        if periodic:
            check_ndarray(pbc, c_bool, 3, "pbc")

        energy = c_double(0.0)
        gradient = np.zeros((natoms, 3), dtype=c_double)
        charges = np.zeros(natoms, dtype=c_double)
        polarizibilities = np.zeros(natoms, dtype=c_double)
        c6_coefficients = np.zeros((natoms, natoms), dtype=c_double)

        if periodic:
            cell_gradient = np.zeros((3, 3), dtype=c_double)
            stress_tensor = np.zeros((3, 3), dtype=c_double)
            args = [
                c_int(natoms),
                numbers.ctypes.data_as(POINTER(c_int)),
                c_double(charge),
                positions.ctypes.data_as(POINTER(c_double)),
                cell.ctypes.data_as(POINTER(c_double)),
                pbc.ctypes.data_as(POINTER(c_bool)),
                output.encode('utf-8'),
                DFTDParameters(**options),
                DFTDOptions(**options),
                energy,
                gradient.ctypes.data_as(POINTER(c_double)),
                None,
                stress_tensor.ctypes.data_as(POINTER(c_double)),
                cell_gradient.ctypes.data_as(POINTER(c_double)),
                polarizibilities.ctypes.data_as(POINTER(c_double)),
                c6_coefficients.ctypes.data_as(POINTER(c_double)),
                charges.ctypes.data_as(POINTER(c_double)),
            ]
            stat = self.library.D4_PBC_calculation(*args)
        else:
            args = [
                c_int(natoms),
                numbers.ctypes.data_as(POINTER(c_int)),
                c_double(charge),
                positions.ctypes.data_as(POINTER(c_double)),
                output.encode('utf-8'),
                DFTDParameters(**options),
                DFTDOptions(**options),
                energy,
                gradient.ctypes.data_as(POINTER(c_double)),
                None,
                polarizibilities.ctypes.data_as(POINTER(c_double)),
                c6_coefficients.ctypes.data_as(POINTER(c_double)),
                charges.ctypes.data_as(POINTER(c_double)),
            ]
            stat = self.library.D4_calculation(*args)

        if stat != 0:
            raise RuntimeError("D4 calculation failed in dftd4.")

        results = {
            'energy': energy.value,
            'gradient': gradient,
            'polarizibilities': polarizibilities,
            'c6_coefficients': c6_coefficients,
            'charges': charges,
        }
        if periodic:
            results['cell gradient'] = cell_gradient
            results['stress tensor'] = stress_tensor
        return results

    def D4Parameters(self, functional: str, lmbd: int = P_MBD_APPROX_ATM) -> dict:
        """load damping parameters from dftd4 shared library"""
        dparam = DFTDParameters()
        stat = self.library.D4_damping_parameters(functional.encode('utf-8'),
                                                  dparam, c_int(lmbd))
        if stat == 0:
            return dparam.to_dict()
        return {}

    # pylint: disable=invalid-name, too-many-arguments, too-many-locals
    def D3Calculation(self, natoms: int, numbers, positions, options: dict,
                      output: str = "-", cell=None, pbc=None) -> dict:
        """wrapper for calling the D3 Calculator from the library."""
        periodic = cell is not None and pbc is not None

        check_ndarray(numbers, c_int, natoms, "numbers")
        check_ndarray(positions, c_double, 3*natoms, "positions")
        if periodic:
            check_ndarray(cell, c_double, 9, "cell")
            check_ndarray(pbc, c_bool, 3, "pbc")

        energy = c_double(0.0)
        gradient = np.zeros((natoms, 3), dtype=c_double)
        polarizibilities = np.zeros(natoms, dtype=c_double)
        c6_coefficients = np.zeros((natoms, natoms), dtype=c_double)

        if periodic:
            cell_gradient = np.zeros((3, 3), dtype=c_double)
            stress_tensor = np.zeros((3, 3), dtype=c_double)
            args = [
                c_int(natoms),
                numbers.ctypes.data_as(POINTER(c_int)),
                positions.ctypes.data_as(POINTER(c_double)),
                cell.ctypes.data_as(POINTER(c_double)),
                pbc.ctypes.data_as(POINTER(c_bool)),
                output.encode('utf-8'),
                DFTDParameters(**options),
                DFTDOptions(**options),
                energy,
                gradient.ctypes.data_as(POINTER(c_double)),
                None,
                stress_tensor.ctypes.data_as(POINTER(c_double)),
                cell_gradient.ctypes.data_as(POINTER(c_double)),
                polarizibilities.ctypes.data_as(POINTER(c_double)),
                c6_coefficients.ctypes.data_as(POINTER(c_double)),
            ]
            stat = self.library.D3_PBC_calculation(*args)
        else:
            args = [
                c_int(natoms),
                numbers.ctypes.data_as(POINTER(c_int)),
                positions.ctypes.data_as(POINTER(c_double)),
                output.encode('utf-8'),
                DFTDParameters(**options),
                DFTDOptions(**options),
                energy,
                gradient.ctypes.data_as(POINTER(c_double)),
                None,
                polarizibilities.ctypes.data_as(POINTER(c_double)),
                c6_coefficients.ctypes.data_as(POINTER(c_double)),
            ]
            stat = self.library.D3_calculation(*args)

        if stat != 0:
            raise RuntimeError("D3 calculation failed in dftd4.")

        results = {
            'energy': energy.value,
            'gradient': gradient,
            'polarizibilities': polarizibilities,
            'c6_coefficients': c6_coefficients,
        }
        if periodic:
            results['cell gradient'] = cell_gradient
            results['stress tensor'] = stress_tensor
        return results
