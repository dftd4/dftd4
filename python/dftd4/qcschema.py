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
"""
QCSchema Support
----------------

Integration with the `QCArchive infrastructure <http://docs.qcarchive.molssi.org>`_.

This module provides a way to translate QCSchema or QCElemental Atomic Input
into a format understandable by the ``dftd4`` API which in turn provides the
calculation results in a QCSchema compatible format.

Supported keywords are

======================== =========== ============================================
 Keyword                  Default     Description
======================== =========== ============================================
 level_hint               None        Dispersion correction level ("d4" or "d4s")
 params_tweaks            None        Optional dict with the damping parameters
 pair_resolved            False       Enable pairwise resolved dispersion energy
 property                 False       Evaluate dispersion related properties
======================== =========== ============================================

The params_tweaks dict contains the damping parameters, at least s8, a1 and a2
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
 ga                       3.0         Charge scaling limiting value
 gc                       2.0         Charge scaling steepness
 wf                       6.0         Coordination number weighting
======================== =========== ============================================

Either method or s8, a1 and a2 must be provided, s9 can be used to overwrite
the ATM scaling if the method is provided in the model.
Disabling the three-body dispersion (s9=0.0) changes the internal selection rules
for damping parameters of a given method and prefers special two-body only
damping parameters if available!

.. note::

    input_data.model.method with a full method name and input_data.keywords["params_tweaks"]
    cannot be provided at the same time. It is an error to provide both options at the
    same time.

Example
-------

>>> from dftd4.qcschema import run_qcschema
>>> import qcelemental as qcel
>>> atomic_input = qcel.models.AtomicInput(
...     molecule = qcel.models.Molecule(
...         symbols = ["O", "H", "H"],
...         geometry = [
...             0.00000000000000,  0.00000000000000, -0.73578586109551,
...             1.44183152868459,  0.00000000000000,  0.36789293054775,
...            -1.44183152868459,  0.00000000000000,  0.36789293054775
...         ],
...     ),
...     driver = "energy",
...     model = {
...         "method": "TPSS-D4",
...     },
...     keywords = {},
... )
...
>>> atomic_result = run_qcschema(atomic_input)
>>> atomic_result.return_result
-0.0002667885779142513
"""

import sys
from typing import Union, overload

import numpy as np

from .interface import DampingParam, DispersionModel
from .library import get_api_version

if sys.version_info < (3, 14):
    try:
        import qcelemental.models.v1 as qcel_v1
    except ModuleNotFoundError:
        import qcelemental.models as qcel_v1
else:
    qcel_v1 = None

try:
    import qcelemental.models.v2 as qcel_v2
except ModuleNotFoundError:
    qcel_v2 = None


if qcel_v1 is None and qcel_v2 is None:
    raise ModuleNotFoundError(
        "The qcelemental package is required for qcschema support. "
        "Please install it with 'pip install qcelemental'."
    )


_supported_drivers = [
    "energy",
    "gradient",
]

_available_levels = [
    "d4",
    "d4s",
]

_clean_dashlevel = str.maketrans("", "", "()")


if qcel_v1 is not None:

    @overload
    def run_qcschema(
        input_data: Union[dict, "qcel_v1.AtomicInput"],
    ) -> "qcel_v1.AtomicResult": ...


if qcel_v2 is not None:

    @overload
    def run_qcschema(
        input_data: Union[dict, "qcel_v2.AtomicInput"],
    ) -> "qcel_v2.AtomicResult": ...


def run_qcschema(input_data):
    """Perform disperson correction based on an atomic inputmodel"""

    if qcel_v2 is not None and isinstance(input_data, qcel_v2.AtomicInput):
        atomic_input = input_data
    elif qcel_v1 is not None and isinstance(input_data, qcel_v1.AtomicInput):
        atomic_input = input_data
    elif qcel_v2 is not None and input_data.get("specification"):
        atomic_input = qcel_v2.AtomicInput(**input_data)
    elif qcel_v1 is not None:
        atomic_input = qcel_v1.AtomicInput(**input_data)
    else:
        raise ValueError(
            "Input data is not a valid QCSchema AtomicInput for either v1 or v2."
        )

    schema_version = atomic_input.schema_version
    if schema_version == 1:
        ret_data = atomic_input.dict()
        input_keywords = atomic_input.keywords
        input_method = atomic_input.model.method
        input_driver = atomic_input.driver
    elif schema_version == 2:
        ret_data = {
            "input_data": atomic_input,
            "extras": {},
            "molecule": atomic_input.molecule,
        }
        input_keywords = atomic_input.specification.keywords
        input_method = atomic_input.specification.model.method
        input_driver = atomic_input.specification.driver
    else:
        raise ValueError(
            f"Unsupported QCSchema version: {schema_version}. Only v1 and v2 are supported."
        )

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
    _level = input_keywords.get("level_hint", "d4")
    if _level.lower() not in _available_levels:
        error = dict(
            error_type="input error",
            error_message="Level '{}' is invalid for this dispersion correction".format(
                _level
            ),
        )
        if schema_version == 1:
            ret_data.update(
                provenance=provenance,
                success=success,
                properties=properties,
                return_result=return_result,
                error=qcel_v1.ComputeError(**error),
            )
            return qcel_v1.AtomicResult(**ret_data)
        elif schema_version == 2:
            return qcel_v2.FailedOperation(
                input_data=atomic_input, error=qcel_v2.ComputeError(**error)
            )

    # Check if the method is provided and strip the “dashlevel” from the method
    _method = input_method.split("-")
    if _method[-1].lower().translate(_clean_dashlevel) == _level.lower():
        _method.pop()
    _method = "-".join(_method)
    if len(_method) == 0:
        _method = None

    # Obtain the parameters for the damping function
    _input_param = input_keywords.get("params_tweaks", {"method": _method})
    if _level.lower() == "d4s":
        _model_param = {
            key: _input_param.pop(key, default)
            for key, default in (
                ("ga", 3.0),
                ("gc", 2.0),
            )
        }
    else:
        _model_param = {
            key: _input_param.pop(key, default)
            for key, default in (
                ("ga", 3.0),
                ("gc", 2.0),
                ("wf", 6.0),
            )
        }

    try:
        param = DampingParam(**_input_param)

        disp = DispersionModel(
            atomic_input.molecule.atomic_numbers[atomic_input.molecule.real],
            atomic_input.molecule.geometry[atomic_input.molecule.real],
            atomic_input.molecule.molecular_charge,
            model=_level,
            **_model_param,
        )

        driver = input_driver
        res = disp.get_dispersion(
            param=param,
            grad=driver == "gradient",
        )
        if input_keywords.get("property", False):
            res.update(**disp.get_properties())
        extras = {"dftd4": res}

        if driver == "gradient":
            if all(atomic_input.molecule.real):
                fullgrad = res.get("gradient")
            else:
                ireal = np.argwhere(atomic_input.molecule.real).reshape((-1))
                fullgrad = np.zeros_like(atomic_input.molecule.geometry)
                fullgrad[ireal, :] = res.get("gradient")

        properties.update(return_energy=res.get("energy"))

        if input_keywords.get("pair_resolved", False):
            res = disp.get_pairwise_dispersion(param=param)
            extras["dftd4"].update(res)

        success = driver in _supported_drivers
        if driver == "energy":
            return_result = properties["return_energy"]
        elif driver == "gradient":
            return_result = fullgrad
        else:
            ret_data.update(
                error=dict(
                    error_type="input error",
                    error_message="Calculation succeeded but invalid driver request provided",
                ),
            )

        ret_data["extras"].update(extras)

    except (RuntimeError, TypeError) as e:
        ret_data.update(
            error=dict(error_type="input error", error_message=str(e)),
        )

    ret_data.update(
        provenance=provenance,
        success=success,
        properties=properties,
        return_result=return_result,
    )

    if schema_version == 1:
        return qcel_v1.AtomicResult(**ret_data)

    if "error" in ret_data:
        return qcel_v2.FailedOperation(
            input_data=atomic_input, error=ret_data["error"]
        )
    return qcel_v2.AtomicResult(**ret_data)
