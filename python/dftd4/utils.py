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
"""Utilities for working with dftd4-data."""

from math import exp
from dftd4.data import def2ecp_nuclear_charges, chemical_hardness


def extrapolate_c8_coeff(c6_coeff: float, iat: int, jat: int) -> float:
    """C6 -> C8 extrapolation"""
    from dftd4.data import sqrt_z_r4_over_r2 as r4r2
    return 3.0 * r4r2[iat] * r4r2[jat] * c6_coeff


def extrapolate_c10_coeff(c6_coeff: float, iat: int, jat: int) -> float:
    """C8 -> C10 extrapolation"""
    c8_coeff = extrapolate_c8_coeff(c6_coeff, iat, jat)
    return 49.0 / 40.0 * c8_coeff * c8_coeff / c6_coeff


def extrapolate_c_n_coeff(c6_coeff: float, iat: int, jat: int, n_order: int)\
        -> float:
    """extrapolation of higher C_n"""
    if n_order == 6:
        return c6_coeff
    if n_order == 8:
        return extrapolate_c8_coeff(c6_coeff, iat, jat)
    if n_order == 10:
        return extrapolate_c10_coeff(c6_coeff, iat, jat)
    if n_order % 2 == 0 and n_order > 10:
        cfrac = extrapolate_c_n_coeff(c6_coeff, iat, jat, n_order-2) / \
                extrapolate_c_n_coeff(c6_coeff, iat, jat, n_order-4)
        return extrapolate_c_n_coeff(c6_coeff, iat, jat, n_order-6) * cfrac**3
    raise ValueError("This is no valid C_n-coefficient")


def cn_gaussian_weight(weighting_factor: float, cni: float, cnref: float) -> float:
    """coordination number based gaussian weight"""
    return exp(-weighting_factor*(cni-cnref)**2)


def gompertz_function(scale: float, gamma: float, zref: float, zmod: float) \
                      -> float:
    """bare scaling function used in dftd4"""
    return exp(scale * (1 - exp(gamma * (1 - zref/zmod))))


def zeta_scaling(g_a: float, g_c: float, charge: float, ref_charge: float,
                 iat: int) -> float:
    """wrapper for scaling function to account for constants and limiting case"""
    zeff = def2ecp_nuclear_charges[iat]
    gamma = g_c * chemical_hardness[iat]
    if zeff+charge > 0.0:
        return gompertz_function(g_a, gamma, zeff+ref_charge, zeff+charge)
    return exp(g_a)
