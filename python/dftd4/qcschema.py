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
"""Integration with the `QCArchive infrastructure <http://docs.qcarchive.molssi.org>`_.

This module provides a way to translate QCSchema or QCElemental Atomic Input
into a format understandable by the ``dftd4`` API which in turn provides the
calculation results in a QCSchema compatible format.

Supported keywords are

======================== =========== ============================================
 Keyword                  Default     Description
======================== =========== ============================================
 level_hint               None        Dispersion correction level (allowed: "d4")
 param_tweaks             None        Optional dict with the damping parameters
======================== =========== ============================================

The param_tweaks dict contains the damping parameters, at least s8, a1 and a2
must be provided

======================== =========== ============================================
 Tweakable parameter      Default     Description
======================== =========== ============================================
 s6                       1.0         Scaling of the dipole-dipole dispersion
 s8                       None        Scaling of the dipole-quadrupole dispersion
 s9                       1.0         Scaling of the three-body dispersion energy
 a1                       None        Scaling of the critical radii
 a2                       None        Offset of the critical radii
 alp                      16.0        Exponent of the zero damping (ATM only)
======================== =========== ============================================

Either method or s8, a1 and a2 must be provided, s9 can be used to overwrite
the ATM scaling if the method is provided in the model.
Disabling the three-body dispersion (s9=0.0) changes the internal selection rules
for damping parameters of a given method and prefers special two-body only
damping parameters if available!
"""

from typing import Union
from .interface import DispersionModel, DampingParam
from .libdftd4 import get_api_version
import qcelemental as qcel


_supported_drivers = [
    "energy",
    "gradient",
]

_available_levels = [
    "d4",
]


def run_qcschema(
    input_data: Union[dict, qcel.models.AtomicInput]
) -> qcel.models.AtomicResult:
    """Perform disperson correction based on an atomic inputmodel"""

    if not isinstance(input_data, qcel.models.AtomicInput):
        atomic_input = qcel.models.AtomicInput(**input_data)
    else:
        atomic_input = input_data
    ret_data = atomic_input.dict()

    provenance = {
        "creator": "dftd4",
        "version": get_api_version(),
        "routine": "dftd4.qcschema.run_qcschema",
    }
    success = False
    return_result = 0.0
    properties = {}

    # Since it is a level hint we a forgiving if it is not present,
    # we are much less forgiving if the wrong level is hinted here.
    _level = atomic_input.keywords.get("level_hint", "d4")
    if _level.lower() not in _available_levels:
        ret_data.update(
            provenance=provenance,
            success=success,
            properties=properties,
            return_result=return_result,
            error=qcel.models.ComputeError(
                error_type="input error",
                error_message="Level '{}' is invalid for this dispersion correction".format(_level),
            ),
        )
        return qcel.models.AtomicResult(**ret_data)

    # Check if the method is provided and strip the “dashlevel” from the method
    _method = atomic_input.model.method
    if len(_method) == 0:
        _method = None
    else:
        _method = _method.split("-")
        if _method[-1].lower() == _level.lower():
            _method.pop()
        _method = "-".join(_method)

    # Obtain the parameters for the damping function
    _input_param = atomic_input.keywords.get("param_tweaks", {})

    try:
        param = DampingParam(
            method=_method,
            **_input_param,
        )

        disp = DispersionModel(
            atomic_input.molecule.atomic_numbers,
            atomic_input.molecule.geometry,
            atomic_input.molecule.molecular_charge,
        )

        res = disp.get_dispersion(
            param=param,
            grad=atomic_input.driver == "gradient",
        )

        properties.update(return_energy=res.get_energy())

        success = atomic_input.driver in _supported_drivers
        if atomic_input.driver == "energy":
            return_result = properties["return_energy"]
        elif atomic_input.driver == "gradient":
            return_result = res.get_gradient()
        else:
            ret_data.update(
                error=qcel.models.ComputeError(
                    error_type="input error",
                    error_message="Calculation succeeded but invalid driver request provided",
                ),
            )

    except RuntimeError as e:
        ret_data.update(
            error=qcel.models.ComputeError(
                error_type="input error", error_message=str(e)
            ),
        ),

    ret_data.update(
        provenance=provenance,
        success=success,
        properties=properties,
        return_result=return_result,
    )

    return qcel.models.AtomicResult(**ret_data)
