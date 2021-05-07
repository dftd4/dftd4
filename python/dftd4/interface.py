# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.
"""Wrapper around the C-API of the dftd4 shared library."""

from typing import Optional
import numpy as np


from .libdftd4 import (
    ffi as _ffi,
    lib as _lib,
    new_error,
    new_structure,
    new_d4_model,
    custom_d4_model,
    new_rational_damping,
    load_rational_damping,
    handle_error,
)


class Structure:
    """
    .. Molecular structure data

    Represents a wrapped structure object in ``dftd4``.
    The molecular structure data object has a fixed number of atoms
    and immutable atomic identifiers
    """

    _mol = _ffi.NULL

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        charge: Optional[float] = None,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
    ):
        """Create new molecular structure data"""
        if positions.size % 3 != 0:
            raise ValueError("Expected tripels of cartesian coordinates")

        if 3 * numbers.size != positions.size:
            raise ValueError("Dimension missmatch between numbers and positions")

        self._natoms = len(numbers)
        _numbers = np.ascontiguousarray(numbers, dtype="i4")
        _positions = np.ascontiguousarray(positions, dtype=float)

        _charge = _ref("double", charge)

        if lattice is not None:
            if lattice.size != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype="float")
        else:
            _lattice = None

        if periodic is not None:
            if periodic.size != 3:
                raise ValueError("Invalid periodicity provided")
            _periodic = np.ascontiguousarray(periodic, dtype="bool")
        else:
            _periodic = None

        self._mol = new_structure(
            self._natoms,
            _cast("int*", _numbers),
            _cast("double*", _positions),
            _charge,
            _cast("double*", _lattice),
            _cast("bool*", _periodic),
        )

    def __len__(self):
        return self._natoms

    def update(
        self,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
    ) -> None:
        """Update coordinates and lattice parameters, both provided in
        atomic units (Bohr).
        The lattice update is optional also for periodic structures.

        Generally, only the cartesian coordinates and the lattice parameters
        can be updated, every other modification, regarding total charge,
        total spin, boundary condition, atomic types or number of atoms
        requires the complete reconstruction of the object.
        """

        if 3 * len(self) != positions.size:
            raise ValueError("Dimension missmatch for positions")
        _positions = np.ascontiguousarray(positions, dtype="float")

        if lattice is not None:
            if lattice.size != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype="float")
        else:
            _lattice = None

        _error = new_error()
        _lib.dftd4_update_structure(
            _error,
            self._mol,
            _cast("double*", _positions),
            _cast("double*", _lattice),
        )

        handle_error(_error)


class DampingParam:
    """Damping parameters for the dispersion correction"""

    _param = _ffi.NULL

    def __init__(self, *, method=None, **kwargs):
        """Create new damping parameter from method name or explicit data"""

        if method is not None:
            _method = _ffi.new("char[]", method.encode())
            self._param = load_rational_damping(
                _method,
                kwargs.get("s9", 1.0) > 0.0,
            )
        else:
            try:
                self._param = new_rational_damping(
                    kwargs.get("s6", 1.0),
                    kwargs["s8"],
                    kwargs.get("s9", 1.0),
                    kwargs["a1"],
                    kwargs["a2"],
                    kwargs.get("alp", 16.0),
                )
            except KeyError as e:
                raise RuntimeError("Constructor requires argument for " + str(e))


class DispersionModel(Structure):
    """
    .. Dispersion model

    Representation of a dispersion model to evaluate C6 coefficients.
    The model is coupled to the molecular structure it has been created
    from and cannot be transfered to another molecular structure without
    recreating it.
    """

    _disp = _ffi.NULL

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        charge: Optional[float] = None,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
        **kwargs,
    ):
        """Create new dispersion model"""

        Structure.__init__(self, numbers, positions, charge, lattice, periodic)

        if "ga" in kwargs or "gc" in kwargs or "wf" in kwargs:
            self._disp = custom_d4_model(
                self._mol,
                kwargs.get("ga", 3.0),
                kwargs.get("gc", 2.0),
                kwargs.get("wf", 6.0),
            )
        else:
            self._disp = new_d4_model(self._mol)

    def get_dispersion(self, param: DampingParam, grad: bool) -> dict:
        """Perform actual evaluation of the dispersion correction"""

        _error = new_error()
        _energy = _ffi.new("double *")
        if grad:
            _gradient = np.zeros((len(self), 3))
            _sigma = np.zeros((3, 3))
        else:
            _gradient = None
            _sigma = None

        _lib.dftd4_get_dispersion(
            _error,
            self._mol,
            self._disp,
            param._param,
            _energy,
            _cast("double*", _gradient),
            _cast("double*", _sigma),
        )

        handle_error(_error)

        results = dict(energy=_energy[0])
        if _gradient is not None:
            results.update(gradient=_gradient)
        if _sigma is not None:
            results.update(virial=_sigma)
        return results

    def get_properties(self) -> dict:
        """Evaluate dispersion related properties"""

        _error = new_error()
        _c6 = np.zeros((len(self), len(self)))
        _cn = np.zeros((len(self)))
        _charges = np.zeros((len(self)))
        _alpha = np.zeros((len(self)))

        _lib.dftd4_get_properties(
            _error,
            self._mol,
            self._disp,
            _cast("double*", _cn),
            _cast("double*", _charges),
            _cast("double*", _c6),
            _cast("double*", _alpha),
        )

        handle_error(_error)

        return {
            "coordination numbers": _cn,
            "partial charges": _charges,
            "c6 coefficients": _c6,
            "polarizibilities": _alpha,
        }

    def get_pairwise_dispersion(self, param: DampingParam) -> dict:
        """Evaluate pairwise representation of the dispersion energy"""

        _error = new_error()
        _pair_disp2 = np.zeros((len(self), len(self)))
        _pair_disp3 = np.zeros((len(self), len(self)))

        _lib.dftd4_get_pairwise_dispersion(
            _error,
            self._mol,
            self._disp,
            param._param,
            _cast("double*", _pair_disp2),
            _cast("double*", _pair_disp3),
        )

        handle_error(_error)

        return {
            "additive pairwise energy": _pair_disp2,
            "non-additive pairwise energy": _pair_disp3,
        }


def _cast(ctype, array):
    """Cast a numpy array to a FFI pointer"""
    return _ffi.NULL if array is None else _ffi.cast(ctype, array.ctypes.data)


def _ref(ctype, value):
    """Create a reference to a value"""
    if value is None:
        return _ffi.NULL
    ref = _ffi.new(ctype + "*")
    ref[0] = value
    return ref
