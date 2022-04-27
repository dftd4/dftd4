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
"""Thin wrapper around the CFFI extension module.

This module mainly acts as a guard for importing the libdftd4 extension and
also provides some FFI based wappers for memory handling.

To check for API compatibility use the provided wrapper around the API version
getter.

Example
-------
>>> from dftd4.library import get_api_version
>>> get_api_version()
'3.4.0'
"""

import functools

try:
    from ._libdftd4 import ffi, lib
except ImportError:
    raise ImportError("dftd4 C extension unimportable, cannot use C-API")


def get_api_version() -> str:
    """Return the current API version from dftd4.
    For easy usage in C the API version is provided as

    10000 * major + 100 * minor + patch

    For Python we want something that looks like a semantic version again.
    """
    api_version = lib.dftd4_get_version()
    return "{}.{}.{}".format(
        api_version // 10000,
        api_version % 10000 // 100,
        api_version % 100,
    )


def _delete_error(error) -> None:
    """Delete a dftd4 error handler object"""
    ptr = ffi.new("dftd4_error *")
    ptr[0] = error
    lib.dftd4_delete_error(ptr)


def new_error():
    """Create new dftd4 error handler object"""
    return ffi.gc(lib.dftd4_new_error(), _delete_error)


def error_check(func):
    """Handle errors for library functions that require an error handle"""

    @functools.wraps(func)
    def handle_error(*args, **kwargs):
        """Run function and than compare context"""
        _err = new_error()
        value = func(_err, *args, **kwargs)
        if lib.dftd4_check_error(_err):
            _message = ffi.new("char[]", 512)
            lib.dftd4_get_error(_err, _message, ffi.NULL)
            raise RuntimeError(ffi.string(_message).decode())
        return value

    return handle_error


def _delete_structure(mol) -> None:
    """Delete molecular structure data"""
    ptr = ffi.new("dftd4_structure *")
    ptr[0] = mol
    lib.dftd4_delete_structure(ptr)


def new_structure(natoms, numbers, positions, charge, lattice, periodic):
    """Create new molecular structure data"""
    return ffi.gc(
        error_check(lib.dftd4_new_structure)(
            natoms,
            numbers,
            positions,
            charge,
            lattice,
            periodic,
        ),
        _delete_structure,
    )


def _delete_model(error) -> None:
    """Delete a dftd4 dispersion model object"""
    ptr = ffi.new("dftd4_model *")
    ptr[0] = error
    lib.dftd4_delete_model(ptr)


def new_d4_model(mol):
    """Create new dftd4 dispersion model object"""
    return ffi.gc(error_check(lib.dftd4_new_d4_model)(mol), _delete_model)


def custom_d4_model(mol, ga, gc, wf):
    """Create new dftd4 dispersion model object"""
    return ffi.gc(
        error_check(lib.dftd4_custom_d4_model)(mol, ga, gc, wf), _delete_model
    )


def _delete_param(error) -> None:
    """Delete a dftd4 damping parameter object"""
    ptr = ffi.new("dftd4_param *")
    ptr[0] = error
    lib.dftd4_delete_param(ptr)


def new_rational_damping(s6, s8, s9, a1, a2, alp):
    """Create new dftd4 damping parameter object"""
    return ffi.gc(
        error_check(lib.dftd4_new_rational_damping)(s6, s8, s9, a1, a2, alp),
        _delete_param,
    )


def load_rational_damping(method, mbd):
    """Create new dftd4 damping parameter object by loading it from the database"""
    return ffi.gc(
        error_check(lib.dftd4_load_rational_damping)(method, mbd),
        _delete_param,
    )


update_structure = error_check(lib.dftd4_update_structure)
get_dispersion = error_check(lib.dftd4_get_dispersion)
get_pairwise_dispersion = error_check(lib.dftd4_get_pairwise_dispersion)
get_properties = error_check(lib.dftd4_get_properties)


def _ref(ctype, value):
    """Create a reference to a value"""
    if value is None:
        return ffi.NULL
    ref = ffi.new(ctype + "*")
    ref[0] = value
    return ref
