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
Tests for the ASE interface to dftd4.
"""

def test_d4_molecular():
    """Short test for DFT-D4 on a subset of the G2."""
    from pytest import approx
    from ase.collections import g2
    from dftd4.calculators import D4_model
    thr = 1.0e-7
    subset = {'C4H4NH': (-0.17254765346459802, 0.0016481180998958495),
              'CH3SCH3': (-0.118396636155385, 0.0013823525532589362),
              'SiH2_s3B1d': (-0.018605813985880824, 0.00015462713987174245),
              'CH3SH': (-0.05722345063877332, 0.0005448462666493717),
              'CH3CO': (-0.05123732582276016, 0.0007742597318805521),
              'C3H4_D2d': (-0.07030306476081483, 0.0009590568552901135),
              'ClF3': (-0.04213170121261114, 0.0007279615294069727),
              'H2CO': (-0.01994226769935472, 0.0002148654202479407),
              'CH3COOH': (-0.08745304310824802, 0.00114445197261909),
              'HCF3': (-0.03887439301151592, 0.00029142286765719783),
              'CH3S': (-0.045072896650802385, 0.0003779814622251629),
              'C3H6_D3h': (-0.09235218055900413, 0.0006339294824490826),
              'CS2': (-0.06617104839166844, 0.001825513384064356),
              'SiH2_s1A1d': (-0.018788257002315355, 0.0003824353060351632),
              'C4H4S': (-0.19935202085994883, 0.0018416124706281756)}

    calc = D4_model(s6=1.0, s8=1.26829475, a1=0.39907098, a2=5.03951304, s9=1.0)
    for name, data in subset.items():
        atoms = g2[name]
        energy, gnorm = data
        atoms.set_calculator(calc)
        assert atoms.get_potential_energy() == approx(energy, thr)
        assert abs(atoms.get_forces().flatten()).mean() == approx(gnorm, thr)
